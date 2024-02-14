import os,re,glob
import numpy as np
import pandas as pd
import argparse
import torch
from transformers import AutoTokenizer, EsmForSequenceClassification


def data_gener(model,tag,outdir,indir):
    ## output directory
    try:
        os.makedirs(outdir)
    except:
        print('Output directory existed: %s.' % (outdir))
        
    ## input directory
    if indir==None:
        file_list = glob.glob('/HOME/data/01.CRISPR/00.data/dataset_0/train'+'/*.faa')
    else:
        file_list = glob.glob(indir+'/*.faa')
    
    ## prepare model
    md = EsmForSequenceClassification.from_pretrained(model, num_labels=4)
    md = md.eval().cuda()
    
    ## get last_hidden_cls
    df = pd.DataFrame(columns=['Set','ProID','Label',] + [ 'c'+str(i) for i in range(0,1280)])
    tokenizer = AutoTokenizer.from_pretrained(model)
    with torch.no_grad():
        for f in file_list:
            label = os.path.basename(f).split('.')[0]
            print('Processing %s ...' % label)
            lines = open(f, 'r').readlines()
            for i in range(len(lines)):
                if lines[i][0] == '>':
                    name = lines[i].strip()[1:]
                    seq = lines[i+1].strip()
                    input = tokenizer(seq, padding=True, truncation=True, max_length=1560, return_tensors='pt', add_special_tokens=True)
                    input = input.to('cuda')
                    output = md(output_hidden_states=True,return_dict=True,**input)
                    emb = output.hidden_states[-1][0][0,:].cpu().tolist()
                    df.loc[len(df.index)] = [tag, name, label] + emb

    ## output last_hidden_cls
    df.to_csv(os.path.join(outdir,tag+'.ft_emb.tab'),sep='\t',index=False)
    torch.cuda.empty_cache()

def main():

    parser = argparse.ArgumentParser(description='data generation')
    parser.add_argument('--model', type=str, default='/HOME/data/01.CRISPR/01.seCasDiscovery/ESM2finetune/multi_class_proka/gamma_0/esm2_t33_650M_UR50D/2023-10-13_09-17-17/epoch_2', help='Fine-tuned model path.')
    parser.add_argument('--tag', type=str, default='train', help="['train','val','test',...]")
    parser.add_argument('--outdir', type=str, help='output directiory')
    parser.add_argument('--indir', type=str, default=None, help='input data directory. default None.')
    args = parser.parse_args()
    
    data_gener(args.model,args.tag,args.outdir,args.indir)
    torch.cuda.empty_cache()

if __name__ == '__main__':

    main()