#!/bin/bash

PROJECT_DIR=$1
TEMP_DIR="${PROJECT_DIR}/runs/${2}_motifs"
CLUSTERING_NAME="liver"

for TISSUE in $PROJECT_DIR/runs/*/; do
	T_NAME=$(basename $TISSUE)   
    for CT in $TISSUE*/; do
	CT_NAME=$(basename $CT)
        cp $CT/motif_discovery_pipeline/3_evaluation/motifs.meme "$TEMP_DIR/${T_NAME}_${CT_NAME}.meme"
    done
done
./getSimilarMotifs.sh $CLUSTERING_NAME $TEMP_DIR/*

#rm -r $TEMP_DIR
