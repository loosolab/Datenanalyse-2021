#!/bin/bash

## Script to cluster all motifs of all runs

# Directory where the whole project lies: 
# It is important that this directory contains the subfolder "runs"
PROJECT_DIR=$1
# Temporary Directory to store MEME files in, Input $2 is the prefix the directory should have
TEMP_DIR="${PROJECT_DIR}/runs/${2}_motifs"
# TODO; File prefix for output file: TODO check if cluster or file? 
CLUSTERING_NAME="liver"

# Temporaray directory to which all the meme files with the motifs will be stored
mkdir $TEMP_DIR

# get all Motifs in meme files of the runs an write them to the temporary directory
for TISSUE in $PROJECT_DIR/runs/*/; do
	T_NAME=$(basename $TISSUE)   
    for CT in $TISSUE*/; do
	    CT_NAME=$(basename $CT)
        if test -f "$CT/motif_discovery_pipeline/3_evaluation/motifs.meme"; then
            cp $CT/motif_discovery_pipeline/3_evaluation/motifs.meme "$TEMP_DIR/${T_NAME}_${CT_NAME}.meme"
        fi
    done
done

# Get path to this script
SPATH=$(dirname $0)
# execute clustering
${SPATH}/utils/getSimilarMotifs.sh $PROJECT_DIR $CLUSTERING_NAME $TEMP_DIR/*

# remove temporary file
rm -r $TEMP_DIR
