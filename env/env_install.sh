# conda create -n chooser python=3.9
# conda activate chooser
conda install pytorch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1 pytorch-cuda=11.7 -c pytorch -c nvidia -y
conda install -c conda-forge tensorboard -y
pip install git+https://github.com/facebookresearch/esm.git
conda install pyg -c pyg -y
conda install -c anaconda scikit-learn -y
conda install -c anaconda scipy -y
pip install git+https://github.com/aqlaboratory/openfold.git
pip install git+https://github.com/huggingface/transformers
conda install -c anaconda ipykernel -y
conda install -c conda-forge einops -y
pip install dm-tree
conda install -c anaconda cudnn -y
conda install -c conda-forge biopython -y
conda install -c conda-forge modelcif -y
conda install -c conda-forge ml-collections -y
conda install -c conda-forge omegaconf -y
conda install -c bioconda bioawk -y
conda install -c conda-forge biotite -y
conda install -c anaconda pandas -y
conda install -c conda-forge matplotlib -y
conda install -c anaconda seaborn -y
conda install -c bioconda samtools -y
conda install -c bioconda seqtk -y
conda install -c bioconda seqkit -y
conda install -c bioconda mafft -y