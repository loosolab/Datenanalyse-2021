#!/bin/bash
# check if conda evironment is active
if  [ ! "$CONDA_DEFAULT_ENV" == "TOBIAS_ENV" ]; then
    echo "Activating TOBIAS_ENV..."
    source /opt/miniconda/bin/activate TOBIAS_ENV
fi

MOTIF_FILES="${@:2}"
NAME=$1
# put all motifs in same File
dt=$(date '+%d-%m-%Y_%H-%M')
joined="${NAME}_motifs_all_${dt}.meme"
TOBIAS FormatMotifs --input $MOTIF_FILES --task join --output $joined --format meme

# apply TOBIAS Clustering
TOBIAS ClusterMotifs --motifs $joined  --clust_method "complete" --dist_method "pcc" --threshold  0.3 --type pdf -o "${NAME}_Cluster"

conda deactivate
