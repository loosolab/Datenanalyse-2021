# Work Package 5 - Motif Discovery

## Introduction

## Setup and Installation
To use this work package it is important that you have three further tools installed:
* [TOBIAS](https://github.com/loosolab/TOBIAS)
* [motif discovery pipeline](https://github.com/loosolab/motif-discovery-pipeline)
* [GSEA](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html)

Please follow the installation instructions of those projects and make sure you have the according conda environments installed. It is neccecary that the TOBIAS environment is called "TOBIAS_ENV". 
Furthermore, you need an environment for snakemake wich should also be called "snakemake". To use the GSEA tool you have to follow the instructions of the "Running GSEA from the Command Line" part.

Additionally 
To install the WP5-Package please ... # TODO
--> TODO: environment with Plotly <br>
To use the GSEA you have to use one of the supported [data formats](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29). With the WP5-Package the *.gmt format can be generated automatically.


## Usage
This work package consists of several scripts (found in the *scripts* directory) which can also be concatenated to a pipeline. For more detailed information please refer to the [WP5-Wiki](https://github.com/loosolab/Datenanalyse-2021/wiki/WP5).

The single steps of the pipeline are:
1. [create_folders.sh](#1-create-folders)
2. [generate configs](#2-generate-configs)

    a) generate_configs.sh
        
    b) generate_configs_with_annotation.sh
    
3. [run_pipeline.py](#3-run-pipeline)
4. [check_logs.sh](#4-check-logs) (optional)
5. [renameAllMotifs.sh](#5-rename-motifs)
6. [cluster_all.sh](#6-cluster-all)
7. [Eval_Motif_similarity.py](#7-evaluate-motif-similarity)

Steps 1-5 are contain the actual motif-discovery-pipeline run, aswell as the needed pre- and postprocessing.
Steps 6 and 7 are for the analysis of similarities between the found motifs in different pipeline runs. 
In case of the given ATAC-Seq data the motif discovery pipelin is ran for each celltype of each tissue given.

Before Starting the pipeline it is important to adjust the config file (TODO).


### 1. Create folders
### 2. Generate configs
### 3. Run pipeline
### 4. Check logs
### 5. Rename motifs
### 6. Cluster all
### 7. Evaluate motif similarity

## Example
TODO
