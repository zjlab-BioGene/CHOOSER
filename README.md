# CHOOSER: Discovering CRISPR-Cas system with self-processing pre-crRNA capability by foundation models

## Introduction

This repository hosts the custom code of **CHOOSER** (**C**as **HO**mlog **O**bserving and **SE**lf-processing sc**R**eening), a novel, effective and unified AI framework for alignment-free discovery of novel CRISPR-Cas systems with self-processing precursor CRISPR RNA (pre-crRNA) capability utilizing protein foundation models. 

![Schematic diagram of the CHOOSER framework for identifying and functional screening of CRISPR-Cas systems with self-processing pre-crRNA capability](Figure_1-2.png)

## Model weights & example data

*Note: before running the notebooks on Colab, please ensure that you have to save/add the following folders into your Google Drive ("MyDrive" folder):

- [model](https://drive.google.com/drive/folders/1y4WKwsoBsqBb_R2Cdj0cwYiLIPnBXj01?usp=sharing)

- [inputs](https://drive.google.com/drive/folders/18GGlIEWYtJVTn2oBXMqbghyYQCLKelLg?usp=sharing)

## Custom code for Colab Notebook

### Step.1 Cas homolog discovery

Firstly, install and run our bioinformatic pipeline [CRISPRCasMiner](https://github.com/zjlab-BioGene/CRISPRCasMiner) to fetch the suspected proteins around the CRISPR arrays [colab_notebook](https://colab.research.google.com/drive/1PYo_vFefUnPWgFLQ5q3Oxu2pTtx9BvzY?usp=sharing):

```
## Run cctyper and prodigal
cctyper example/input_test.fna output/01_cctyper \
    --db data \
    --prodigal meta \
    --keep_tmp

## Run CRISPRCasMiner
python ccminer/ccminer.py example/input_test.fna output/02_ccminer \
    --cctyper_path output/01_cctyper \
    --database_name my_project_name \
    --name my_sample_name \
    --db data \
    --prodigal meta \
    --span 10 \
    --keep_tmp
```

Secondly, 

01_CasDiscovery [colab_notebook](https://colab.research.google.com/drive/1oxa1YrmgCe5ok7GwWCuHwGoZ1M_Otikr?usp=sharing)

`Input`: suspicious proteins in .fasta/.faa format.

`Output`: predicted tags of the proteins, including `cas9`, `cas12`, `cas13` and `other`.



### Step.2 Pre-crRNA self-processing functional screening of Cas12 homologs 
02_Cas12_SelfProcessing [colab_notebook](https://colab.research.google.com/drive/1D5_Qffq-EUZYQk_tTKMftCv9wvxSh2Kz?usp=sharing)

`Input`: Cas12 proteins in .fasta/.faa format.
`Output`: predicted wether the Cas12 candidates would self-process their pre-crRNAs.

## Contacts

qiliu@tongji.edu.cn, liwh@zhejianglab.com, liwh@tongji.edu.cn