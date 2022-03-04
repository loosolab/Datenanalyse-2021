#!/bin/bash
# Path to project directory
DIR=$1

for TISSUE in $DIR/*/; do
    TIS_NAME=$(basename $TISSUE)
    for CT in $TISSUE*/; do
	    CT_NAME=$(basename $CT)
        PATTERN="${TIS_NAME}_${CT_NAME}"
        echo $CT
        sed -E -i "s/motif_([[:digit:]]+)([[:blank:]]+)motif_([[:digit:]]+)/${PATTERN}_\1\2${PATTERN}_\3/" $CT/motif_discovery_pipeline/3_evaluation/motifs.meme
    done
done