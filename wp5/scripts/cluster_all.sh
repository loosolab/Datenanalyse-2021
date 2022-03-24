#!/bin/bash

## Script to cluster all motifs of all runs

# get script path
SPATH=$(dirname $0)
# read in config
CONF="${SPATH}/../global_vars.cnf"
while read LINE; do declare "$LINE"; done < $CONF

# Temporary Directory to store MEME files in, Input $2 is the prefix the directory should have
TEMP_DIR="${PROJECT_DIR}/runs/tmp_motifs_dir"
# TODO; File prefix for output file: TODO check if cluster or file? 
CLUSTERING_NAME=$1

# Create temporaray directory to which all the MEME files with the motifs will be stored
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

# execute clustering
${SPATH}/utils/getSimilarMotifs.sh $PROJECT_DIR $CLUSTERING_NAME $TEMP_DIR/*

# remove temporary file
rm -r $TEMP_DIR
