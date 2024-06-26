# isoFarm

A Snakemake pipeline for RNA isoforms quantification with Salmon, RSEM and DEXseq.

## Installation

Clone the repo on your machine
```
git clone git@github.com:molinerisLab/isoFarm.git isoFarm
```
Create conda environment and activate it before running the workflow
```
conda env create -f local/env/environment.yml -n isoFarm
conda activate isoFarm
```

## Usage

Move to `dataset/v1` and set project specific configuration by modifying the `config.yaml` (understand [Fragment Library Type](https://salmon.readthedocs.io/en/latest/library_type.html#fraglibtype)).

Then run
```
snakemake -p -j N_CORES all
```

To visualize the isoforms of a GENE_NAME run
```
snakemake -p -j N_CORES Salmon_plot/GENE_NAME_transcripts.bar.pdf
```
