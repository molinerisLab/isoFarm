# isoFarm

A Snakemake pipeline for RNA isoforms quantification with Salmon, RSEM and DEXseq.


## Dependencies

All software dependencies required to run this pipeline are specified in the env.yaml file located local/envs.

To set up the environment, you only need to have either `Conda` or `Mamba` installed on your system.

### IPAFinder and rMATS

If you plan to use differential [IPAFinder](https://github.com/ZhaozzReal/IPAFinder) or [rMATS](https://github.com/Xinglab/rmats-turbo) functionalities, you must manually clone the corresponding external repository in a shared directory `/path/to/external-repo`, as it is not included directly in this pipeline to avoid redundancy of external data.

Once cloned, the location of the external repository should be specified in your `config.yaml` file using the appropriate key:

```
EXTERNAL_REPO_DIR: "/path/to/external-repo"
```

Make sure the path points to a directory that contains the required scripts and data used by differential IPA and rMATS modules.


## Installation

Clone the repo on your machine
```
git clone git@github.com:molinerisLab/isoFarm.git isoFarm
cd isoFarm
```
Create conda environment with all the dependencied and activate it
```
conda env create -f local/env/environment.yaml -n isoFarm_Env
conda activate isoFarm_Env
```

## Usage

Manually create a `metadata.txt` file (tab separated) having the sample names in the first column (mandatory name: `sample`), and the sample annotations in the other columns.

Move to `dataset/v1` and set project specific configuration by modifying the `config.yaml`, then run one of the following according to your needs.

```
snakemake -p -j N all_salmon
snakemake -p -j N all_rsem
snakemake -p -j N all_ipafinder
```

To identify the top expressed isoform for condition and all the expressed isoforms in each condition, set the `METADATA_GROUP_COL` in the `config.yaml` and run:

```
snakemake -p -j N rsem.isoforms.results.gene_symbol.avg_condition.top_exp_isoform.gz rsem.isoforms.results.gene_symbol.avg_condition.expressed.gz
```

