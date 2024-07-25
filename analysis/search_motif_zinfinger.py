import os,sys,re
import argparse
import shutil
import pandas as pd
import gzip
from Bio import SeqIO

motif = r'C.{2}C' # Cas13 HEPN domain, conserved R(N)xxxxH motif

def search_fasta(infa,out):
    try:
        os.mkdir(out)
        print('Ouput fold %s created!' % out)
    except:
        print('Ouput fold %s exists!' % out)
    fasta = gzip.open(fasta,'rt') if '.gz' in infa else infa
    outab = open(out+'/zinfinger_mtf.tab','w')
    outab.write('\t'.join(['protein','plen','mtf_num','position','motifs','sequence'])+'\n')
    i = 0 # number of proteins
    j = 0 # number of proteins with motif matches
    for record in list(SeqIO.parse(fasta,'fasta')):
        i = i+1
        seq = str(record.seq)
        length = str(len(seq))
        pid = record.id
        mtf_match = list(re.finditer(motif,seq))
        mtf_num = len(mtf_match)
        if mtf_num > 0:
            j = j+1
            mtfs = []
            poses = []
            for mtf in mtf_match:
                mtf_str = mtf.group(0)
                pos = str(mtf.span()).replace(' ','')
                mtfs.append(mtf_str)
                poses.append(pos)
            outab.write('\t'.join([pid,length,str(mtf_num),';'.join(poses),';'.join(mtfs),seq])+'\n')
    outab.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--protein',type=str,default=None,required=True,help='input protein fasta file')
    parser.add_argument('--out',type=str,default=None,required=True,help='output name')
    args = parser.parse_args()

    search_fasta(args.protein,args.out)