import os, re
import time
import shutil
import ast
import logging
import pandas as pd
import numpy as np
import pyfaidx
from sklearn.utils import shuffle
import torch
from torch.utils.data import DataLoader, Dataset
from tqdm import tqdm
import argparse
from transformers import AutoTokenizer, EsmForSequenceClassification

## Class
class TqdmLoggingHandler(logging.Handler):
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except Exception:
            self.handleError(record)

class SequenceDataset(Dataset):
    def __init__(self, inputs, labels, names):
        self.input_ids = inputs['input_ids']
        self.attention_mask = inputs['attention_mask']
        self.labels = torch.tensor(labels)
        self.names = names
    def __len__(self):
        return len(self.input_ids)
    def __getitem__(self, idx):
        return {'labels': self.labels[idx], 'input_ids': self.input_ids[idx], 'attention_mask': self.attention_mask[idx], 'ids': idx}
    def get_num_samples_per_class(self):
        return torch.bincount(self.labels).tolist()

## Constants
DATASET_TRAINING_KEYS = ['labels', 'input_ids', 'attention_mask']
all_labels = ['cas9','cas12','cas13','noncas']
dataset_random_seed = None

## Functions
def parse_head_mask(value):
    try:
        print(torch.tensor(ast.literal_eval(value)))
        return torch.tensor(ast.literal_eval(value))
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid head_mask value: {value}")

def create_dataset(tokenizer, fasta_dir, max_seq_len, label_to_id_fn, random_seed):
    labels, sequences, names = read_fasta(fasta_dir)
    if random_seed is not None:
        labels, sequences, names = shuffle(labels, sequences, names, random_state=random_seed) # type: ignore
    inputs = tokenizer(sequences, padding=True, truncation=True, max_length=max_seq_len, return_tensors='pt', add_special_tokens=True)
    label_ids = label_to_id_fn(labels)
    return SequenceDataset(inputs, label_ids, names)

def read_fasta(fasta_dir):
    labels = []
    names = []
    sequences = []
    for fasta_file in os.listdir(fasta_dir):
        if not fasta_file.endswith(('.faa', '.fasta')):
            continue
        label = fasta_file.split('.')[0]
        fasta = pyfaidx.Fasta(os.path.join(fasta_dir, fasta_file), rebuild=False)
        for record in fasta:
            labels.append(label)
            seq = str(record)
            seq = re.sub(r"[\n\*]", '', seq)
            seq = re.sub(r"[UZOB]", "X", seq)
            sequences.append(seq)
            names.append(record.name)
    print(f"Read {len(labels)} sequences from {fasta_dir}, sequences: {len(sequences)}, names: {len(names)} from fasta_dir: {fasta_dir}")
    time.sleep(1) # avoid multi process issues
    return labels, sequences, names

def get_dataloader(tokenizer, label_to_id_fn, random_seed, args):
    eval_dataset = create_dataset(tokenizer, args.eval_dataset_dir, args.max_seq_len, label_to_id_fn, random_seed)
    eval_dataloader = DataLoader(eval_dataset, batch_size=args.eval_batch_size)
    return eval_dataloader

def get_label_to_id_fn(all_labels):
    def label_to_id_fn(labels):
        return [all_labels.index(label) if label in all_labels else 0 for label in labels]
    return label_to_id_fn

def get_parameters():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', '-i', type=str, default=None, help='Input suspicious proteins in fasta.')
    ap.add_argument('--output', '-o', type=str, default='pred_result.csv', help='Output prediction table.')
    ap.add_argument('--model', '-m', type=str, default='./models/CasDiscovery', help='Pre-trained model path.')
    ap.add_argument('--max_len', type=int, default=1560, help='Max-sequence-length (aa).')
    ap.add_argument('--batch_size', type=int, default=1, help='Batch size.')
    ap.add_argument('--seed', type=int, default=42, help='Random seed.')
    args = ap.parse_args()
    if args.input == None or args.model == None:
        raise Exception('Error: input suspicious protein.faa and pre-tained model are required!')
    args.max_seq_len = args.max_len
    args.eval_batch_size = args.batch_size
    
    args.eval_dataset_dir = './protein'
    try:
        shutil.rmtree(args.eval_dataset_dir)
        os.makedirs(args.eval_dataset_dir)
    except:
        os.makedirs(args.eval_dataset_dir)
    os.system('cp %s %s/suspicious.faa' % (args.input,args.eval_dataset_dir ))
    
    return args

def run_predict(args):
    label_to_id_fn = get_label_to_id_fn(all_labels)
    tokenizer = AutoTokenizer.from_pretrained(args.model)
    model = EsmForSequenceClassification.from_pretrained(args.model, num_labels=len(all_labels)).cuda().eval()
    eval_dataloader = get_dataloader(tokenizer, label_to_id_fn, dataset_random_seed, args)

    merged_ids = torch.tensor([], dtype=torch.int)
    merged_predicted_label_ids = torch.tensor([], dtype=torch.int)
    merged_logits = torch.tensor([], dtype=torch.float)
    
    with torch.no_grad():
        for batch in eval_dataloader:
            inputs = {k: v for k, v in batch.items() if k in DATASET_TRAINING_KEYS}
            label_ids = batch.get("labels")
            ids = batch.get("ids")
            inputs['labels'] = inputs['labels'].cuda()
            inputs['input_ids'] = inputs['input_ids'].cuda()
            inputs['attention_mask'] = inputs['attention_mask'].cuda()
            outputs = model(**inputs)
            logits = outputs.get("logits")
            predicted_probs = torch.softmax(logits, dim=1)
            predicted_label_ids = torch.argmax(logits, dim=1)

            merged_ids = torch.cat((merged_ids, ids.cpu()), dim=0)
            merged_predicted_label_ids = torch.cat((merged_predicted_label_ids, predicted_label_ids.cpu()), dim=0)
            merged_logits = torch.cat((merged_logits, predicted_probs.cpu()), dim=0)

            merged_predicted_labels = [all_labels[label_id] for label_id in merged_predicted_label_ids]
            merged_names = [eval_dataloader.dataset.names[id] for id in merged_ids]
            df = pd.DataFrame({
                "name": merged_names,
                "predicted_label": merged_predicted_labels,
                })
            for i, label in enumerate(all_labels):
                df[f"prob: {label}"] = [f"{round(prob[i] * 100, 2)}%" for prob in merged_logits.numpy()]
                df.to_csv(args.output, sep='\t',index=False)

    print('Prediction finished!\nResults output to: %s' % args.output)
    print('------------------------------')
    print('Total protein: %d' % df.shape[0])
    print('Cas9: %d' % df[df.predicted_label=='cas9'].shape[0])
    print('Cas12: %d' % df[df.predicted_label=='cas12'].shape[0])
    print('Cas13: %d' % df[df.predicted_label=='cas13'].shape[0])

## Run
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(TqdmLoggingHandler())

## main   
def main():
    args = get_parameters()
    run_predict(args)
    
## -------------------------------- Run -------------------------------- ##
if __name__ == '__main__':
    main()  

# Environment: +pyfaidx
## pip install pyfaidx

# Run:

## export CUDA_VISIBLE_DEVICES=0
# python CasDiscovery.py -m YOUR_MODEL_PATH -i YOUR_INPUT.suspicious.faa -o OUTPUT_PREDICTION.suspicious.csv