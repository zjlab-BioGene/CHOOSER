import os, sys, re
import pandas as pd
import numpy as np
import argparse
import datetime
import logging
import shutil
import subprocess
import multiprocess as mp
import glob
import tqdm

logging.basicConfig(
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level='INFO'
)

def get_fname(file):
    fname = 'Unknown'
    if file.endswith('pdb'):
        fname = re.sub(r'\.pdb$', '', file)
    elif file.endswith('pdb.gz'):
        fname = re.sub(r'\.pdb\.gz$', '', file)
    elif file.endswith('cif'):
        fname = re.sub(r'\.cif$', '', file)
    fname = re.sub(':', '|', fname)
    return fname

def _tuple2str(tup):
    s, e = tup
    return f"{s},{e}"

def _slice_seq(tup, seq):
    s, e = tup
    return seq[s-1:e]

def _concat_to_regions_worker(alignpos, protein, cutoff=20):
    # 连续对齐区间（exactly aligned residues）
    region_pos = []
    start = 1
    end = 1
    # find first aligned position
    for i in range(0, len(alignpos)):
        if alignpos[i] > 0:
            start = i + 1
            end = i + 1
            break
    for i in range(start - 1, len(alignpos)):
        if alignpos[i] > 0:
            if alignpos[i - 1] > 0:
                end = i + 1
            else:
                start = i + 1
        else:
            if alignpos[i - 1] > 0:
                end = i
                region_pos.append((start, end))
            else:
                continue
    if (start, end) not in region_pos:
        region_pos.append((start, end))

    region_pos_info = ';'.join([_tuple2str(x) for x in region_pos])
    region1 = [_slice_seq(x, protein) for x in region_pos]

    # 合并 gap 小于等于 cutoff 的区间
    region_pos2 = []
    start2, end2 = region_pos[0]
    start2 = int(start2); end2 = int(end2)
    for i in range(1, len(region_pos)):
        last_s, last_e = region_pos[i - 1]
        this_s, this_e = region_pos[i]
        if this_s - last_e <= cutoff:
            end2 = this_e
        else:
            region_pos2.append((start2, end2))
            start2 = this_s
            end2 = this_e
    if (start2, end2) not in region_pos2:
        region_pos2.append((start2, end2))

    region2 = [_slice_seq(x, protein) for x in region_pos2]
    region2_seq = ''.join(region2)
    region_pos2_info = ';'.join([_tuple2str(x) for x in region_pos2])
    return region_pos_info, region_pos, region1, region_pos2_info, region_pos2, region2_seq, region2

def _parse_one_alignment(args):
    aln_path, cutoff = args
    aln_base = os.path.basename(aln_path)
    q_name, db_name = re.sub(r'\.aln$', '', aln_base).split('@')[0:2]
    remove_align = 0

    with open(aln_path, 'r') as handle:
        lines = handle.readlines()

    # TM-scores（第 14/15 行）
    tmscore_1 = float(lines[14].split()[1])
    tmscore_2 = float(lines[15].split()[1])
    if (tmscore_1 < 0.6) and (tmscore_2 < 0.6):
        remove_align = 1

    # 长度（行中 L=）
    db_len = int(re.sub('L=', '', lines[14].split()[7].strip(',')))
    q_len  = int(re.sub('L=', '', lines[15].split()[7].strip(',')))

    # 序列（20/21 行）
    atag = re.sub('\n', '', lines[20])
    pseq = re.sub('\n', '', lines[21])
    con = list(re.sub(' ', 'X', atag))
    tag = list(np.zeros(len(con), dtype=int))
    for i in range(len(con)):
        if con[i] != 'X':
            tag[i] += 1

    protein = re.sub('-', '', pseq)

    # 将对齐位置映射回原 query（不含 '-'）
    alignpos = list(np.zeros(q_len, dtype=int))
    region_seq = ''
    pos = 0
    q_seq = list(pseq)
    for i in range(len(q_seq)):
        if i > len(tag) - 1:
            break
        aa = q_seq[i]
        if aa != '-':
            pos += 1
            if tag[i] > 0:
                region_seq += aa
                alignpos[pos - 1] += 1

    region_pos_info, region_pos, region1, region_pos2_info, region_pos2, region2_seq, region2 = _concat_to_regions_worker(
        alignpos, protein, cutoff=cutoff
    )

    # per-position counts
    region_counts = list(np.zeros(len(protein), dtype=int))
    for (s, e) in region_pos2:
        for i in range(s - 1, e):
            region_counts[i] += 1

    # region_report rows
    region_report_rows = []
    for i, (s, e) in enumerate(region_pos2):
        pos_str = f"{s},{e}"
        seg = region2[i]
        region_report_rows.append([
            q_name, db_name, q_len, db_len, tmscore_2, tmscore_1, i + 1, pos_str, len(seg), seg
        ])

    # summary row
    summary_row = [
        q_name, db_name, q_len, db_len,
        tmscore_2, tmscore_1,
        region_pos_info, len(region_seq), region_seq,
        region_pos2_info, len(region2_seq), region2_seq
    ]
    return {
        'remove': remove_align,
        'aln_path': aln_path,
        'q_name': q_name,
        'protein': protein,
        'summary_row': summary_row,
        'region_report_rows': region_report_rows,
        'region_counts': region_counts,
    }

def _append_tsv(path, rows, columns):
    """Append rows to a TSV file; create with header if not exists."""
    if not rows:
        return
    df = pd.DataFrame(rows, columns=columns)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    exists = os.path.exists(path)
    df.to_csv(path, sep='\t', index=False, mode='a', header=not exists)

def _process_alignment(args):
    # 运行 USalign → 解析 → 删除 .aln → 按需删除空 .err
    db_root, in_root, out_dir, q, db, env_threads, cutoff, keep_err = args
    q_name = get_fname(q); db_name = get_fname(db)
    aln_path = os.path.join(out_dir, f"{q_name}@{db_name}.aln")
    err_path = os.path.join(out_dir, f"{q_name}@{db_name}.err")

    env = os.environ.copy()
    for var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[var] = str(env_threads)

    # 1) Run USalign（允许非 0 返回码，不中断流水线）
    with open(aln_path, 'w') as out_f, open(err_path, 'a') as err_f:
        subprocess.run(
            ['USalign', os.path.join(db_root, db), os.path.join(in_root, q)],
            stdout=out_f, stderr=err_f, env=env, check=False
        )

    # 2) 如果 .aln 缺失或空，返回一个最小记录
    if not os.path.exists(aln_path) or os.path.getsize(aln_path) == 0:
        parsed = {
            'remove': 1,
            'aln_path': aln_path,
            'q_name': get_fname(q),
            'protein': '',
            'summary_row': None,
            'region_report_rows': [],
            'region_counts': [],
        }
    else:
        # 3) 解析
        try:
            parsed = _parse_one_alignment((aln_path, cutoff))
        except Exception:
            parsed = {
                'remove': 1,
                'aln_path': aln_path,
                'q_name': get_fname(q),
                'protein': '',
                'summary_row': None,
                'region_report_rows': [],
                'region_counts': [],
            }

    # 4) 无论成功与否，都删除 .aln
    try:
        if os.path.exists(aln_path):
            os.remove(aln_path)
    except Exception:
        pass

    # 5) 若指定 --keep_err，则删除**空**的 .err（仅保留非空文件）
    try:
        if keep_err and os.path.exists(err_path) and os.path.getsize(err_path) == 0:
            os.remove(err_path)
    except Exception:
        pass

    return parsed

# ------------------ CLI ------------------
default_db = '/data/wenhui.li/00.database/HMMdatabase/StructureAlign/RuvC_proteins'
ap = argparse.ArgumentParser()
ap.add_argument('--input','-i',type=str,default=None,required=True,help='Input query pdb directory. required.')
ap.add_argument('--db','-d',type=str,default=default_db,required=True,help='Input database directory. required.')
ap.add_argument('--out','-o',type=str,default='./',required=True,help='Output directory. required.')
ap.add_argument('--cpu','-t',type=int,default=4,required=False,help='Processes. default 4.')
ap.add_argument('--cutoff','-l',type=int,default=20,required=False,help='Maximum insertion length to merge regions. default 20.')
ap.add_argument('--resume', action='store_true', help='Resume: if .aln exists and non-empty, parse instead of re-running USalign.')
ap.add_argument('--chunksize', type=int, default=32, required=False, help='Chunksize for imap_unordered (IPC efficiency).')
ap.add_argument('--per_query_subdir', action='store_true', help='Put outputs in subdirs per query to reduce FS contention.')
ap.add_argument('--env_threads', type=int, default=1, help='Set OMP/MKL/OPENBLAS threads per USalign process.')
ap.add_argument('--keep_err', action='store_true', help='Keep only non-empty .err files (empty .err will be deleted).')
ap.add_argument('--ignore_lst', type=str, default=None, help='Ignore list file. If specified, the queries in the file will be ignored.')
ap.add_argument('--suffix', type=str, default='pdb', help='Suffix of ignore list files. default pdb.')
ap.add_argument('--start_method', type=str, choices=['spawn','fork','forkserver'], default='spawn',
                help='Multiprocessing start method. Default spawn (safer for OpenMP/BLAS).')

# ------------------ Batch Pipeline ------------------

def get_ignore_lst(ignore_lst, suffix):
    if ignore_lst is None:
        return []
    with open(ignore_lst, 'r') as f:
        return [ '.'.join([x.strip(), suffix]) for x in f.readlines() ]

class BatchUSalign(object):
    def __init__(self, args):
        self.input = args.input
        self.suffix = args.suffix
        self.db = args.db
        self.out = args.out
        self.ignore_lst = get_ignore_lst(args.ignore_lst, self.suffix)
        self.queries = [x for x in os.listdir(self.input) if x.endswith(('.pdb', '.pdb.gz', '.cif')) and x not in self.ignore_lst]
        self.db_files = [x for x in os.listdir(self.db) if x.endswith(('.pdb', '.pdb.gz', '.cif'))]
        self.db_num = len(self.db_files)
        self.qdb = [(q, db) for q in self.queries for db in self.db_files]

        self.out_align = os.path.join(self.out, 'usalign')
        # self.out_report = os.path.join(self.out, 'report')

        self.cpu = args.cpu
        self.cutoff = args.cutoff
        self.resume = args.resume
        self.per_query_subdir = args.per_query_subdir
        self.env_threads = args.env_threads
        self.keep_err = args.keep_err
        self.chunksize = args.chunksize

        # 汇总容器
        self.proteins = {}     # q_name -> protein sequence
        self.region = {}       # q_name -> per-position counts
        self.pair = {}

        # 输出 DataFrame（内存中保留，以便最后兜底写出）
        self.summary = pd.DataFrame(columns=[
            'Query','Database','q_len','db_len','qTMscore','dbTMscore',
            'Exactly_region','Exactly_len','Exactly_seq',
            'regionPos','regionLen','regionSeq'
        ])
        self.region_report = pd.DataFrame(columns=[
            'Query','Database','q_len','db_len','qTMscore','dbTMscore',
            'Rank','regionRange','Rlen','Segment'
        ])

        logging.basicConfig(
            format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            level='INFO'
        )

        # 初始化输出目录
        if os.path.exists(self.out):
            if not self.resume:
                shutil.rmtree(self.out)
                os.makedirs(self.out, exist_ok=True)
            else:
                os.makedirs(self.out, exist_ok=True)
        else:
            os.makedirs(self.out, exist_ok=True)
        os.makedirs(self.out_align, exist_ok=True)
        # os.makedirs(self.out_report, exist_ok=True)

        # 流式写出目标
        self.summary_path = os.path.join(self.out, 'REPORT_USalign.txt')
        if not self.resume:
            if os.path.exists(self.summary_path):
                os.remove(self.summary_path)

        # per-query 完成计数（用于释放内存）；键统一使用 get_fname(q)
        self.q_total = {get_fname(q): self.db_num for q in self.queries}
        self.q_done = {get_fname(q): 0 for q in self.queries}

    def run_batch(self):
        self.run_search()
        self.write_summary()

    def run_search(self):
        logging.info('Running USalign against database...')
        tasks = []
        for q, db in self.qdb:
            q_name = get_fname(q)
            out_dir = self.out_align if not self.per_query_subdir else os.path.join(self.out_align, q_name)
            os.makedirs(out_dir, exist_ok=True)
            aln_path = os.path.join(out_dir, f"{q_name}@{get_fname(db)}.aln")

            # resume 策略：若 .aln 已存在且非空，则只解析并删除
            if self.resume and os.path.exists(aln_path) and os.path.getsize(aln_path) > 0:
                tasks.append((self.db, self.input, out_dir, q, db, self.env_threads, self.cutoff, self.keep_err))
                continue

            tasks.append((self.db, self.input, out_dir, q, db, self.env_threads, self.cutoff, self.keep_err))

        if not tasks:
            logging.info('Nothing to do (all alignments present).')
            return

        total_tasks = len(tasks)
        # --- Adaptive chunksize: avoid starving workers when total_tasks is small ---
        eff_chunksize = self.chunksize if (self.chunksize and self.chunksize > 0) else 1
        # If user-specified chunksize is too large relative to total tasks and workers,
        # reduce it so that we have at least ~4 chunks per worker.
        target_chunks = max(4 * self.cpu, 1)
        if total_tasks < eff_chunksize * self.cpu:
            eff_chunksize = max(1, total_tasks // target_chunks)
        if eff_chunksize <= 0:
            eff_chunksize = 1
        logging.info(f"Scheduler: total_tasks={total_tasks}, workers={self.cpu}, chunksize={eff_chunksize} (requested={self.chunksize})")

        ctx = mp.get_context(self.start_method if hasattr(self, 'start_method') else 'spawn')
        with ctx.Pool(processes=self.cpu, maxtasksperchild=64) as pool:
            with tqdm.tqdm(total=total_tasks, desc=f"USalign→parse→stream ({self.cpu} procs, chunk={eff_chunksize})", smoothing=0.1) as pbar:
                for result in pool.imap_unordered(_process_alignment, tasks, chunksize=eff_chunksize):
                    # 若 worker 无返回，直接推进进度
                    if not result:
                        pbar.update(1)
                        continue

                    # 统一使用 qn（get_fname 之后的名称）计数；若不存在则初始化
                    qn = result.get('q_name')
                    if qn is not None:
                        if qn not in self.q_done:
                            # 兼容：某些情况下 qn 可能未在初始化字典中（例如外部变换了文件名）
                            self.q_total[qn] = self.db_num
                            self.q_done[qn] = 0
                        self.q_done[qn] += 1
                        if self.q_done[qn] == self.q_total[qn]:
                            # 该 query 的所有 db 对齐结束，释放内存（若存在）
                            if qn in self.region:
                                del self.region[qn]
                            if qn in self.proteins:
                                del self.proteins[qn]

                    # 被过滤/解析失败则不做下游聚合，只更新进度
                    if result.get('summary_row') is None:
                        pbar.update(1)
                        continue

                    protein = result['protein']
                    if qn not in self.proteins:
                        self.proteins[qn] = protein

                    # 聚合逐位计数
                    rc = result['region_counts']
                    if qn not in self.region:
                        self.region[qn] = list(np.zeros(len(self.proteins[qn]), dtype=int))
                    if len(self.region[qn]) != len(rc):
                        logging.warning(f"Length mismatch for {qn}: stored {len(self.region[qn])} vs new {len(rc)}; using min length.")
                    L = min(len(self.region[qn]), len(rc))
                    for i in range(L):
                        self.region[qn][i] += rc[i]

                    # 立刻写出 REPORT_USalign.txt（并保留内存副本兜底）
                    _append_tsv(self.summary_path, [result['summary_row']], self.summary.columns.tolist())
                    self.summary.loc[len(self.summary.index)] = result['summary_row']

                    # # 立刻写出 per-query report（并保留内存副本兜底）
                    # if result['region_report_rows']:
                    #     q_report_path = os.path.join(self.out_report, f"{qn}_alignReport.txt")
                    #     _append_tsv(q_report_path, result['region_report_rows'], self.region_report.columns.tolist())
                    #     for row in result['region_report_rows']:
                    #         self.region_report.loc[len(self.region_report.index)] = row

                    pbar.update(1)

    # 兼容保留（未使用，但保留旧 API）
    def parse_batch_usalign(self):
        logging.info("parse_batch_usalign is bypassed (now parsing on-the-fly during run_search).")
        return

    def write_summary(self):
        # 兜底：若仍有缓冲未落盘（极少），再写一次
        if len(self.summary) > 0:
            _append_tsv(self.summary_path, self.summary.values.tolist(), self.summary.columns.tolist())
        logging.info('Streaming write completed (SUMMARY_REPORT disabled): %s', self.summary_path)

def main():
    args = ap.parse_args()
    master = BatchUSalign(args)
    master.run_batch()

if __name__ == '__main__':
    # 显式设置 start_method，避免在某些 shell/notebook 环境下的兼容问题
    try:
        mp.set_start_method(ap.parse_args().start_method, force=True)
    except RuntimeError:  
        pass  
    main()