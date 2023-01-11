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

Move to `dataset/v1` and set project specific configuration by modifying the `config.yaml` file.

Then run
```
snakemake -p -j N_CORES all
```
