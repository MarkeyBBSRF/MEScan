# MEScan for large scale ME patterns survey
We provide an example dataset from TCGA Ovarian Cancer to demonstrate the pipeline for large scale ME patterns survey using MEScan. 

![pipeline](/mcmc/snakemake/dag.mcmc.png)

Dependencies:
  * R:
    1. MEScan
    2. locfdr
    3. snowfall
  * Snakemake

Execuate the following code to executate the pipeline. 
```
cd mcmc
snakemake -s snake_me/Snakefile -p re/ov/ov.hc_sets.txt
```

# How to generate MEScan input
MEScan uses the input generated by R package `MADGiC`. We provide the source codes and data to create the input file for MEScan.
## Prepare mutation table
The first step convert the MAF table to a

## Create MESCan input
Download the required data and source codes from [here](link_to_madgic). 

The folder contains the following files:
```
# MAGGiC functions
MADGiC.R
# MADGiC input data
exome_36.RData
gene.rep.expr.RData
gene_names.txt
prior.RData
```

The input can be created using `create_mescan_input.R` script providing the mutation table created in the previous step and the desired output folder.
```
Rscript create_mescan_input.R demo/mutationtable.Rdata demo_output
```

It will take some time to run the input preparation pipeline.