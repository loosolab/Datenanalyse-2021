#!/bin/bash

## Script to rename all motifs of the runs with the naming schema Tissue_CellType_Nr
## Caution: Only the meme files containing all the motifs will be changed. The JASPAR file containing single Motifs will not be changed.

# get script path
SPATH=$(dirname $0)
# read in config
CONF="${SPATH}/../global_vars.cnf"
while read LINE; do declare "$LINE"; done < $CONF

# Path to "runs" directory
DIR="${PROJECT_DIR}/runs"
SUCCESS=true

# Iterate over all tissues and their celltypes and change the motif names.
echo "Renaming motifs"
for TISSUE in $DIR/*/; do
    TIS_NAME=$(basename $TISSUE)
    for CT in $TISSUE*/; do
	    CT_NAME=$(basename $CT)
        # check if file exists
        FILE=${CT}/motif_discovery_pipeline/3_evaluation/motifs.meme
        if test -f "$FILE"; then
            PATTERN="${TIS_NAME}-${CT_NAME}"
	    {
            sed -E -i "s/motif_([[:digit:]]+)([[:blank:]]+)motif_([[:digit:]]+)/${PATTERN}-\1\2${PATTERN}-\3/" $FILE
    	    } || {
	    SUCCESS=false    	    
	    }
    			
        fi
        
    done
done

if [ "$SUCCESS" = true ] ; then
	echo "Renaming Successful"
	exit 0
fi
echo "Renaming failed for at least one file"
exit 1
