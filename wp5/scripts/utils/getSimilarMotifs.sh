#!/bin/bash

## Script to cluster all motifs of given files

# check if conda evironment is active
if  [ ! "$CONDA_DEFAULT_ENV" == "TOBIAS_ENV" ]; then
    echo "Activating TOBIAS_ENV..."
    source /opt/miniconda/bin/activate TOBIAS_ENV
fi

# Threshhold for clustering
THRESH=$1
# Directory where to write the output to
DIR=$2
# File prefix for output file:
NAME=$3
# All Motif files to cluster
MOTIF_FILES="${@:4}"

# put all motifs in same File
dt=$(date '+%d-%m-%Y_%H-%M')
joined="${DIR}/${NAME}_motifs_all_${dt}.meme"
TOBIAS FormatMotifs --input $MOTIF_FILES --task join --output $joined --format meme

# apply TOBIAS Clustering
TOBIAS ClusterMotifs --motifs $joined  --clust_method "complete" --dist_method "pcc" --threshold  $THRESH --type pdf -o "${DIR}/${NAME}_Cluster" 

conda deactivate
