#!/bin/bash
# check if conda evironment is active
if  [ ! "$CONDA_DEFAULT_ENV" == "TOBIAS_ENV" ]; then
    echo "Please activate TOBIAS_ENV environment first"
    echo "Abort..."
    exit 1
fi

MOTIF_DIR="${@:1}"

# put all motifs in same File
dt=$(date '+%d-%m-%Y_%H-%M')
joined="joined_motifs_all_${dt}.jaspar"
TOBIAS FormatMotifs --input MOTIF_DIR --task join --output $joined

# apply TOBIAS Clustering
# TODO adjust hyper parameters
TOBIAS ClusterMotifs --motifs $joined --threshold  0.3 --type pdf -o Similar_Motifs_TOBIAS