#!/bin/bash

## Script to rename all motifs of the runs with the naming schema Tissue_CellType_Nr
## Caution: Only the meme files containing all the motifs will be changed. The JASPAR file containing single Motifs will not be changed.

# Path to project directory
DIR=$1

# Iterate over all tissues and their celltypes and change the motif names.
echo "Renaming motifs"
for TISSUE in $DIR/*/; do
    TIS_NAME=$(basename $TISSUE)
    for CT in $TISSUE*/; do
	    CT_NAME=$(basename $CT)
        PATTERN="${TIS_NAME}_${CT_NAME}"
        sed -E -i "s/motif_([[:digit:]]+)([[:blank:]]+)motif_([[:digit:]]+)/${PATTERN}_\1\2${PATTERN}_\3/" $CT/motif_discovery_pipeline/3_evaluation/motifs.meme
    done
done

echo "Renaming done!"