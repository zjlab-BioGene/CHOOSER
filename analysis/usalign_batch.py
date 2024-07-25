import os,sys,re
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

logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level='INFO')

## Parameters
default_db = '/data/wenhui.li/00.database/HMMdatabase/StructureAlign/RuvC_proteins'
ap = argparse.ArgumentParser()
ap.add_argument('--input','-i',type=str,default=None,required=True,help='Input query datatable. required.')
ap.add_argument('--db','-d',type=str,default=default_db,required=True,help='Input database directory. required.')
ap.add_argument('--out','-o',type=str,default='./',required=True,help='Output Directory. required.')
ap.add_argument('--cpu','-t',type=int,default=4,required=False,help='Threads. default 4.')
ap.add_argument('--cutoff','-l',type=int,default=20,required=False,help='Maximum insertion length. default 20.')

## BatchData Process
class BatchUSalign(object):
    def __init__(self, args):
        self.input = args.input
        self.queries = os.listdir(self.input)
        self.db = args.db
        self.db_num = len(os.listdir(self.db))
        self.out = args.out
        self.out_align = os.path.join(self.out,'usalign')
        self.out_report = os.path.join(self.out,'report')
        self.cpu = args.cpu
        self.cutoff = args.cutoff
        self.proteins = {}     
        self.qdb = []
        for q in self.queries:
            for db in os.listdir(self.db):
                self.qdb.append((q,db))
        self.region = {}
        self.pair = {}
                
        self.summary = pd.DataFrame(columns=['Query','Database','q_len','db_len','qTMscore','dbTMscore','Exactly_region','Exactly_len','Exactly_seq','regionPos','regionLen','regionSeq'])
        self.region_report = pd.DataFrame(columns=['Query','Database','q_len','db_len','qTMscore','dbTMscore','Rank','regionRange','Rlen','Segment'])
        self.summary_report = pd.DataFrame(columns=['Query','q_len','Threadhold','regionPos','regionLen','regionSeq'])
                
        logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level='INFO')
        
        try:
            shutil.rmtree(self.out)
        except:
            print("Output dir created: %s" % self.out)
        os.makedirs(self.out)
        os.makedirs(self.out_align)
        os.makedirs(self.out_report)
    
    def run_batch(self):
        self.run_search()
        self.parse_batch_usalign()
        self.write_summary()
    
    def run_search(self):
        logging.info('Running USalign against database...')
        def usalign(qdb):
            def get_fname(file):
                fname = 'Unknown'
                if file.endswith('pdb'):
                    fname = re.sub('\.pdb$','',file)
                elif file.endswith('cif'):
                    fname = re.sub('\.cif$','',file)
                # fname = re.sub('\.','_',fname)
                fname = re.sub(':','|',fname)
                return fname
            query,db = qdb
            db_name = get_fname(db)
            q_name = get_fname(query)
            if (q_name+'@'+db_name not in self.pair.keys()) and (db_name+'@'+q_name not in self.pair.keys()):
                self.pair[q_name+'@'+db_name] = 1
                logging.debug('Running USalign, query: %s, db: %s' % (q_name,db_name))
                with open(os.path.join(self.out,'usalign.log'),'a') as usalign_log:
                    with open(os.path.join(self.out_align, q_name + '@' + db_name+'.aln'),'w') as usalin_out:
                        subprocess.run(['USalign', os.path.join(self.db,db), os.path.join(self.input,query)],
                                    stdout=usalin_out, 
                                    stderr=usalign_log)
        ## Start multiprocess
        pool = mp.Pool(self.cpu)
        # list(tqdm.tqdm(pool.imap(self.usalign, self.qdb), total=len(self.qdb)))
        with tqdm.tqdm(total=len(self.qdb), desc="Processing") as pbar:
            for i in pool.imap(usalign, self.qdb):
                pbar.update(1)
        pool.close()
        pool.join()
    
    def slice_seq(self,tup,seq):
        s,e = tup
        slc = seq[s-1:e]
        return slc  
    
    def tuple2str(self,tup):
        s,e = tup
        return '%d,%d' % (s,e)
    
    def concat_to_regions(self,alignpos,protein,threadshold=0):
        region_pos = []
        start = 1
        end = 1
        for i in range(0,len(alignpos)):
            if alignpos[i] > 0:
                start = i+1
                end = i+1
                break
                
        for i in range(start-1,len(alignpos)):
            if alignpos[i] > threadshold:
                if alignpos[i-1] > threadshold:
                    end = i+1
                else:
                    start = i+1
            else:
                if alignpos[i-1] > threadshold:
                    end = i
                    region_pos.append((start,end))
                else:
                    continue
        if (start,end) not in region_pos:
            region_pos.append((start,end))
        region_pos_info = ';'.join([self.tuple2str(x) for x in region_pos]) ## exactly aligned positions
        region1 = [self.slice_seq(x,protein) for x in region_pos]
        
        ## fetch aligned region
        region_pos2 = []
        start2,end2 = region_pos[0]
        start2 = int(start2)
        end2 = int(end2)
        for i in range(1,len(region_pos)):
            last_s, last_e = region_pos[i-1]
            this_s, this_e = region_pos[i]
            if this_s - last_e <= self.cutoff:
                end2 = this_e
            else:
                region_pos2.append((start2,end2))
                start2 = this_s
                end2 = this_e
                
        if (start2,end2) not in region_pos2:
            region_pos2.append((start2,end2))
        
        region2 = [self.slice_seq(x,protein) for x in region_pos2] ## aligned regions
        region2_seq = ''.join(region2)
        region_pos2_info = ';'.join([self.tuple2str(x) for x in region_pos2]) ## aligned region positions
        
        return region_pos_info,region_pos,region1,region_pos2_info,region_pos2,region2_seq,region2
     
    def stat_pos(self,q_name,protein,region_pos2):
        if q_name not in self.region.keys():
            self.region[q_name] = list(np.zeros(len(protein),dtype=int))
        for (s,e) in region_pos2:
            for i in range(s-1,e):
                self.region[q_name][i] += 1     
    
    def parse_seq(self,q_name,q_len,atag,pseq):
        atag = re.sub('\n','',atag)
        pseq = re.sub('\n','',pseq)
        con = list(re.sub(' ','X',atag))
        tag = list(np.zeros(len(con),dtype=int))
        for i in range(len(con)):
            if con[i] != 'X':
                tag[i] += 1
        # print(tag)
        protein = re.sub('-','',pseq)
        if q_name not in self.proteins.keys():
            self.proteins[q_name] = protein
        q_seq = list(pseq)
        ## fetch aligned residues
        region_seq = '' ## exactly aligned residues
        alignpos = list(np.zeros(q_len,dtype=int))
        pos = 0
        for i in range(len(q_seq)):
            if i > len(tag)-1: 
                break
            aa = q_seq[i]
            if aa != '-':
                pos += 1
                if tag[i] > 0:
                    region_seq += aa
                    alignpos[pos-1] += 1
        
        
        region_pos_info,region_pos,region1,region_pos2_info,region_pos2,region2_seq,region2 = self.concat_to_regions(alignpos,protein,threadshold=0)
        self.stat_pos(q_name,protein,region_pos2)
        
        return region_pos_info,region_pos,region_seq,region1,region_pos2_info,region_pos2,region2_seq,region2
        
    def parse_one_usalign(self,aln):
        q_name,db_name = aln.split('@')[0:2]
        with open(os.path.join(self.out_align,aln),'r') as handle:
            lines = handle.readlines()
            ## TM-scores
            tmscore_1 = float(lines[14].split()[1])
            tmscore_2 = float(lines[15].split()[1])
            ## Sequence length
            db_len = int(re.sub('L=','',lines[14].split()[7].strip(',')))
            q_len = int(re.sub('L=','',lines[15].split()[7].strip(',')))
            ## Parce aligned sequence
            region_pos_info,region_pos,region_seq,region1,region_pos2_info,region_pos2,region2_seq,region2 = self.parse_seq(q_name,q_len,lines[20],lines[21])
            ## save data
            self.summary.loc[len(self.summary.index)] = [q_name,db_name,q_len,db_len,tmscore_2,tmscore_1,region_pos_info,len(region_seq),region_seq,region_pos2_info,len(region2_seq),region2_seq]
            for i in range(len(region_pos2)):
                s,e = region_pos2[i]
                pos = '%d,%d' % (s,e)
                self.region_report.loc[len(self.region_report.index)] = [q_name,db_name,q_len,db_len,tmscore_2,tmscore_1,i+1,pos,len(region2[i]),region2[i]]
    
    def parse_batch_usalign(self):
        for aln in os.listdir(self.out_align):
            self.parse_one_usalign(aln)
        
        for q_name in self.proteins.keys():
            for i in range(0,11):
                thd = int(i*self.db_num/10)
                region_pos_info,region_pos,region1,region_pos2_info,region_pos2,region2_seq,region2 = self.concat_to_regions(self.region[q_name],self.proteins[q_name],threadshold=thd)
                self.summary_report.loc[len(self.summary_report.index)] = [q_name,len(self.proteins[q_name]),i*10,region_pos2_info,len(region2_seq),region2_seq]
                
    def write_summary(self):
        ## write Summary Report
        summary_report = self.summary_report.sort_values(by=['Query','Threadhold'],ascending=[True,True])
        summary_report.to_csv(os.path.join(self.out,'SUMMARY_REPORT_USalign.txt'),sep='\t',index=False)
        ## write REPORT
        summary = self.summary.sort_values(by=['Query','Database'])
        summary.to_csv(os.path.join(self.out,'REPORT_USalign.txt'),sep='\t',index=False)
        ## write each single report
        for qname,df in self.region_report.groupby(by='Query'):
            df = df.sort_values(by=['Database','Rank'],ascending=[True,True])
            df.to_csv(os.path.join(self.out_report,qname+'_alignReport.txt'),sep='\t',index=False)           

## Run
master = BatchUSalign(ap.parse_args())
master.run_batch()
