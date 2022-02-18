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
Running multi_fasta_mask
```
snakemake --cores all --config input={input} 
```

# Options
In the form ```--config <x>=<userx> <y>=<usery>```

Required:

```input```:  Should point to a .fasta formatted file containing multiple sequeces 

Optional

```output_dir```: The Dir that will be created to hold the results of the run. (Default = results)

```AF_THRESHOLD```: The allele frequency threshold for the position to be masked (Default = 0, 0<=x<=1) All varients included 
