#!/bin/bash

## Script to cluster all motifs of given files

# check if conda evironment is active
if  [ ! "$CONDA_DEFAULT_ENV" == "TOBIAS_ENV" ]; then
    echo "Activating TOBIAS_ENV..."
    source /opt/miniconda/bin/activate TOBIAS_ENV
fi

# File prefix for output file: TODO check if cluster or file? 
NAME=$1
# All Motif files to cluster
MOTIF_FILES="${@:2}"

# put all motifs in same File
dt=$(date '+%d-%m-%Y_%H-%M')
joined="${NAME}_motifs_all_${dt}.meme"
TOBIAS FormatMotifs --input $MOTIF_FILES --task join --output $joined --format meme

# apply TOBIAS Clustering
TOBIAS ClusterMotifs --motifs $joined  --clust_method "complete" --dist_method "pcc" --threshold  0.3 --type pdf -o "${NAME}_Cluster"
# TODO delete joinded? 

conda deactivate
