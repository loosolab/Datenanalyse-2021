#!/bin/bash

# input parameters
DIR="/mnt/workspace_stud/allstud/wp5/runs"
FILE_PATH="/mnt/workspace_stud/allstud/wp5/runs"  # put into the evaluation folder

# output file
FILE_NAME=motif_counts.txt
FILE=${FILE_PATH}/$FILE_NAME
if ! [ -f $FILE ]; then
    cd $FILE_PATH
    touch $FILE_NAME
    echo "$FILE_NAME was created."
    echo -e "tissue\tcell_type\tmotif\t" >> $FILE
fi

# count motifs per tissue per cell type and write them into the output file
for TISSUE in $DIR/*; do
    if [ -d $TISSUE ]; then
        TISSUE=$(echo $TISSUE | rev | cut -d'/' -f-1 | rev)
        TISSUE_AVAILABLE=$(cat $FILE | grep -c "${TISSUE}")
        if [ "$TISSUE_AVAILABLE" = 0  ]; then
            for CELL_TYPE in $DIR/$TISSUE/*; do
                if [ -d $CELL_TYPE/motif_discovery_pipeline ]; then
                    CELL_TYPE=$(echo $CELL_TYPE | rev | cut -d'/' -f-1 | rev)
                    COUNTER_MOTIFS=0
                    if [ -d $DIR/$TISSUE/$CELL_TYPE/motif_discovery_pipeline/2_discovery/2_processed_motif/motifs ]; then
                        for MOTIF in $DIR/$TISSUE/$CELL_TYPE/motif_discovery_pipeline/2_discovery/2_processed_motif/motifs/*.png; do
                            let COUNTER_MOTIFS++
                        done
                    fi
                    echo -e "${TISSUE}\t${CELL_TYPE}\t${COUNTER_MOTIFS}\t" >> $FILE
                fi
            done
        else
            echo "$TISSUE is already analyzed."
        fi
    fi
done
