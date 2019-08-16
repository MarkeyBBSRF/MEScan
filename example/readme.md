# MEScan for large scale ME patterns survey
We provide an example dataset from TCGA Ovarian Cancer to demonstrate the pipeline for large scale ME patterns survey using MEScan. The pipeline has several steps as shown in the DAG below.

![pipeline](mcmc/snakemake/dag.mcmc.png)

Dependencies:
  * R:
    1. MEScan
    2. locfdr
    3. snowfall
  * Snakemake



The script below shows how to run MEScan on the example data. The pipeline is managed by the Snakemake workflow language. It reads the rules predefined in `snakemake/Snakefile` file and executes all the required steps to run the MEScan survey. When successes, the high-confidence mutually exclusive gene sets can be found in ` re/ov/ov.hc_sets.genes.txt`.

```
cd mcmc
snakemake -s snakemake/Snakefile -p re/ov/ov.hc_sets.genes.txt
```

# MEScan input
MEScan requires two tables: 1) mutation table and 2) background mutation rate table as inputs for the large scale survey. The inputs are generated using the methods published in the R package **[MADGIC](link_to_madgic) ** from  a given Mutation Annotation Format (MAF) file. We provide a script that take the MAF file as input and output the required input data matrix for MEScan.

#### a. Download MAGDIC related files

The functions and dependencies of MAGIC can be downloaded from [Google Drive](https://drive.google.com/a/g.uky.edu/file/d/15-QjT-8hqWin8gtJQzi5Y-wW3808_UKr/view?usp=sharing). The folder contain the following data and functions. Copy them to the `example/prep_input` folder.

```
MADGiC.R
exome_36.RData
gene.rep.expr.RData
gene_names.txt
prior.RData
```
#### b. Create MEScan input

The `create_mescan_input.R ` takes a mutation annotation file (MAF) as input and creates the input data required by MEScan in `RData` format. 

```
Rscript create_mescan_input.R example.maf output_dir
```



# MEScan output

A list of high-confidence mutually exclusive gene sets will be returned in the following format. The first column is the gene sets delimited by `;`. The second column is the TG score for the gene set.

| Genesets                                        | TG    |
| :---------------------------------------------- | :---- |
| ENSG00000142273;ENSG00000196531;ENSG00000187837 | 10.2  |
| ENSG00000177885;ENSG00000170961;ENSG00000165916 | 9.448 |

