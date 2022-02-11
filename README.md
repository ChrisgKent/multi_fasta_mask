# multi_fasta_mask
A pipeline for the generation of Masked genomes from a multi_fasta file

To install
```
git clone https://github.com/ChrisgKent/multi_fasta_mask
```
Activate / install the conda enviroment
```
cd multi_fasta_mask
conda env create -f multi_fasta_mask.yaml

conda activate multifa_mask
```
Clone required repo from github
```
cd multi_fasta_mask
git clone https://github.com/ChrisgKent/ncov-random-scripts
```
Running multi_fasta_mask
```
snakemake --cores all --config input={input} 
```
