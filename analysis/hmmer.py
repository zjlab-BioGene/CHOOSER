import os,sys,re
import subprocess
import logging
import glob
import tqdm
import multiprocess as mp
import pandas as pd
import argparse
import shutil

ref_hmm_pdir = '/home/wenhui.li/01.data/00.database/Cas/Profiles'

class HMMER():
    
    def __init__(self):
        self._parser = argparse.ArgumentParser()
    
    def parameters(self):
        self._parser.add_argument('--protein',type=str,default=None,required=True,help='Proteins in fasta format. Required.')
        self._parser.add_argument('--out',type=str,default=None,required=True,help='Output directory. Required.')
        self._parser.add_argument('--hmmdb',type=str,default=ref_hmm_pdir,help='Hmm database for protein annotation.')
        self._parser.add_argument('--cpu',type=int,default=1,help='Number of thresholds.')
        self.args = self._parser.parse_args()

        if (self.args.protein == None or self.args.out == None):
            raise Exception('Error: input protein fasta and output directory is required!')

    def hmmsearch(self,hmm): # A single search
        hmm_name = re.sub('\.hmm', '', hmm)   
        logging.debug('Running HMMER against '+hmm_name)
        hmm_out = os.path.join(self.args.out+'hmmer', hmm_name+'.tab')
        with open(self.args.out+'hmmer.log', 'a') as hmmer_log:
            subprocess.run(['hmmsearch', 
                            '--domtblout', hmm_out, 
                            os.path.join(self.args.hmmdb, hmm), 
                            self.args.protein], 
                            stdout=subprocess.DEVNULL, 
                            stderr=hmmer_log)

    ## Step 1, run hmmsearch
    def run_hmm(self):  
        logging.info('Running HMMER against custom profiles')
        # Make dir
        try:
            shutil.rmtree(self.args.out+'hmmer')
        except:
            print("Output dir created: %s" % self.args.out+'hmmer')
        os.makedirs(self.args.out+'hmmer')
        # Start multiprocess
        pool = mp.Pool(self.args.cpu)
        list(tqdm.tqdm(pool.imap(self.hmmsearch, os.listdir(self.args.hmmdb)), total=len(os.listdir(self.args.hmmdb))))
        pool.close()

    ## Step 2, load hmmer reuslts and write out
    def load_and_write_hmm(self):
        logging.debug('Loading HMMER output')  
        # Get files
        hmm_files = glob.glob(os.path.join(self.args.out+'hmmer', '*.tab'))  
        # Parse externally  
        with open(self.args.out+'hmmer.tab', 'w') as hmmer_tab:
            subprocess.run(['grep', '-v', '^#']+hmm_files, stdout=hmmer_tab)
            subprocess.run(['sed', '-i', 's/:/ /', self.args.out+'hmmer.tab'])

        # Load
        hmm_df = pd.read_csv(self.args.out+'hmmer.tab', sep='\s+', header=None,
            usecols=(0,1,3,6,7,
                        8,16,17,18,19,
                        20,21,22),
            names=('Hmm','ORF','tlen','qlen','Eval',
                    'score','hmm_from','hmm_to','ali_from','ali_to',
                    'env_from','env_to','pprop'))

        # Parse HMM names
        hmm_df['Hmm'] = [re.sub('\.tab', '', 
                        re.sub(os.path.join(self.args.out, 'hmmer', ''), '', x)) 
                        for x in hmm_df['Hmm']]
    
        # Add columns
        hmm_df['Gene'] = ['|'.join(x.split('|')[0:-1]) for x in hmm_df['ORF']]
        hmm_df['Sample'] = [x.split('|').pop() for x in hmm_df['ORF']]

        # Coverages of aligments
        def covs(df_sub):
            df_sub['Cov_seq'] = len(set([x for sublst in [list(range(i,j)) 
                for i,j in zip(df_sub['ali_from'], df_sub['ali_to']+1)] 
                for x in sublst])) / df_sub['tlen']
            df_sub['Cov_hmm'] = len(set([x for sublst in [list(range(i,j)) 
                for i,j in zip(df_sub['hmm_from'], df_sub['hmm_to']+1)] 
                for x in sublst])) / df_sub['qlen']
            df_sub = df_sub[['Hmm','ORF','tlen','qlen','Eval','score',
                            'Sample','Gene','Cov_seq','Cov_hmm']]
            df_sub = df_sub.drop_duplicates()
            return df_sub

        hmm_df = hmm_df.groupby(['Hmm','ORF']).apply(covs)
        hmm_df.reset_index(drop=True, inplace=True)
        self.hmm_df = hmm_df.drop_duplicates()
        self.hmm_df.to_csv(self.args.out+'hmmer.tab', sep='\t', index=False)

def main():
    hmmer = HMMER()
    hmmer.parameters()
    hmmer.run_hmm()
    hmmer.load_and_write_hmm()

if __name__ == '__main__':
    sys.exit(main())

