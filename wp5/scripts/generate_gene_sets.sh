#!/bin/bash

# script to generate txt files for the gene set analyzation

# get script path
SPATH=$(dirname $0)
# read in config
CONF="${SPATH}/../test.conf"
while read LINE; do declare "$LINE"; done < $CONF

# input parameters
DIR="${PROJECT_DIR}/runs"
FILE_PATH=$DIR  # same as directory

# generate and fill txt files for the gene set analyzation
for TISSUE in $DIR/*; do
    if [ -d $TISSUE ]; then
        TISSUE=$(echo $TISSUE | rev | cut -d'/' -f-1 | rev) # extract the tissue name
        FILE_PATH_TISSUE="${FILE_PATH}/$TISSUE"
        DT_TISSUE=$(date '+%d/%m/%Y %H:%M:%S')   # start time of the tissue
        echo "Starting ${TISSUE} at ${DT_TISSUE}."
        
        for CELL_TYPE in $FILE_PATH_TISSUE/*; do
            if [ -d $CELL_TYPE ]; then
                CELL_TYPE=$(echo $CELL_TYPE | rev | cut -d'/' -f-1 | rev)    # extract the cell type name
                FILE_PATH_CELL_TYPE="${FILE_PATH_TISSUE}/$CELL_TYPE"
                # generate the file for the distance up to 1000 bp with uniq genes
                FILE_NAME_1K=${TISSUE}_${CELL_TYPE}_1k.txt
                FILE_PATH_SAVE_1K="${FILE_PATH_CELL_TYPE}/annotation"   # path to save the file
                FILE_1K="${FILE_PATH_SAVE_1K}/$FILE_NAME_1K"
                if ! [ -f $FILE_1K ]; then  # only generate file if it isn't already existing
                    cd $FILE_PATH_SAVE_1K
                    touch $FILE_NAME_1K
                    echo "$FILE_NAME_1K was created."
                fi
                # generate the file for the distance up to 2000 bp with unique genes
                FILE_NAME_2K=${TISSUE}_${CELL_TYPE}_2k.txt
                FILE_PATH_SAVE_2K="${FILE_PATH_CELL_TYPE}/annotation"   # path to save the file
                FILE_2K="${FILE_PATH_SAVE_2K}/$FILE_NAME_2K"
                if ! [ -f $FILE_2K ]; then  # only generate file if it isn't already existing
                    cd $FILE_PATH_SAVE_2K
                    touch $FILE_NAME_2K
                    echo "$FILE_NAME_2K was created."
                fi
                # generate the file for the distance up to 1000 bp with all genes
                FILE_NAME_1K_ALL=${TISSUE}_${CELL_TYPE}_1k_all.txt
                FILE_PATH_SAVE_1K_ALL="${FILE_PATH_CELL_TYPE}/annotation"   # path to save the file
                FILE_1K_ALL="${FILE_PATH_SAVE_1K_ALL}/$FILE_NAME_1K_ALL"
                if ! [ -f $FILE_1K_ALL ]; then  # only generate file if it isn't already existing
                    cd $FILE_PATH_SAVE_1K_ALL
                    touch $FILE_NAME_1K_ALL
                    echo "$FILE_NAME_1K_ALL was created."
                fi
                # generate the file for the distance up to 2000 bp with all genes
                FILE_NAME_2K_ALL=${TISSUE}_${CELL_TYPE}_2k_all.txt
                FILE_PATH_SAVE_2K_ALL="${FILE_PATH_CELL_TYPE}/annotation"   # path to save the file
                FILE_2K_ALL="${FILE_PATH_SAVE_2K_ALL}/$FILE_NAME_2K_ALL"
                if ! [ -f $FILE_2K_ALL ]; then  # only generate file if it isn't already existing
                    cd $FILE_PATH_SAVE_2K_ALL
                    touch $FILE_NAME_2K_ALL
                    echo "$FILE_NAME_2K_ALL was created."
                fi
                
                FILE_PATH_ANNOTATION="${FILE_PATH_CELL_TYPE}/motif_discovery_pipeline/4_annotation/uropa"
                
                # check if the cell type is already analyzed
                MOTIF_COUNTER=0  # counter for the to be analyzed motifs
                PUNCTUATION_CHECKER=$(cat $FILE_1K | grep -c "#")   # counts the analyzed motifs
                # counts the to be analyzed motifs
                if [ -d $FILE_PATH_ANNOTATION ]; then
                    for MOTIF in $FILE_PATH_ANNOTATION/*; do
                        MOTIF=$(echo $MOTIF | rev | cut -d'/' -f-1 | rev)    # extract the motif name
                        if [[ $MOTIF != "config" ]];then
                            let MOTIF_COUNTER++
                        fi
                    done
                fi
                # compare the counts
                if [ $MOTIF_COUNTER -gt $PUNCTUATION_CHECKER ]; then
                    DT_CELL_TYPE=$(date '+%d/%m/%Y %H:%M:%S')   # start time of the cell type
                    echo "Starting ${CELL_TYPE} at ${DT_CELL_TYPE}."
                    for MOTIF in $FILE_PATH_ANNOTATION/*; do
                        MOTIF=$(echo $MOTIF | rev | cut -d'/' -f-1 | rev)    # extract the motif name
                        # check if the TF has already been analyzed
                        MOTIF_CHECKER=$(cat $FILE_1K | grep -c $MOTIF)   #  check if the TF is in the file
                        if [ $MOTIF_CHECKER -eq 0 ]; then
                            if [[ $MOTIF != "config" ]]; then
                                DT_MOTIF=$(date '+%d/%m/%Y %H:%M:%S')   # start time of the motif
                                echo "Starting ${MOTIF} at ${DT_MOTIF}."
                                
                                GENE_SET_NAME_1K="#${MOTIF}_1k" # column 1 -> 1k
                                GENE_SET_NAME_2K="#${MOTIF}_2k" # column 1 -> 2k
                                GENES_1K=()                     # columns after 2 -> 1k
                                GENES_2K=()                     # columns after 2 -> 2k
                                    
                                cd $FILE_PATH_ANNOTATION/$MOTIF
                                    
                                # read *allhits.txt
                                # check if the feature = gene $7, check the distance $12, add the gene name $18 to array
                                ALLHITS=($(find . | grep allhits.txt | sed 's,./,,g'))
                                GENES_1K=($(cat $ALLHITS | awk -F'\t' '{if ($7 == "gene" && $12 != "NA" && $12 <= 1000 && $18 != "NA") print $18;}'))
                                GENES_2K=($(cat $ALLHITS | awk -F'\t' '{if ($7 == "gene" && $12 != "NA" && $12 <= 2000 && $18 != "NA") print $18}'))
                                
                                # only sort genes -> no unique genes
                                SORTED_GENES_1K=($(echo "${GENES_1K[@]}" | tr ' ' '\n' | sort | tr '\n' ' '))
                                SORTED_GENES_2K=($(echo "${GENES_1K[@]}" | tr ' ' '\n' | sort | tr '\n' ' '))
                                
                                # only choose the unique genes
                                SORTED_UNIQUE_GENES_1K=($(echo "${GENES_1K[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
                                SORTED_UNIQUE_GENES_2K=($(echo "${GENES_2K[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
                                                                
                                # counter for indices in arrays (SORTED_UNIQUE_GENES_1K)
                                COUNTER_ARR_1K=0
                                # fill the *1k.txt file with data
                                echo -e "${GENE_SET_NAME_1K}" >> $FILE_1K   # insert the gene set name to the file (1k)
                                # append the genes to the assciated gene set name (1k)
                                for elem in "${SORTED_UNIQUE_GENES_1K[@]}"; do
                                    echo -e "${SORTED_UNIQUE_GENES_1K[COUNTER_ARR_1K]}" >> $FILE_1K
                                    let COUNTER_ARR_1K++
                                done
                                echo -e "" >> $FILE_1K    # insert a new line at the end of the motif data

                                # counter for indices in arrays (SORTED_UNIQUE_GENES_2K)
                                COUNTER_ARR_2K=0
                                # fill the *2k.txt file with data
                                echo -e "${GENE_SET_NAME_2K}" >> $FILE_2K   # insert the gene set name to the file (2k)
                                # append the genes to the assciated gene set name (2k)
                                for elem in "${SORTED_UNIQUE_GENES_2K[@]}"; do
                                    echo -e "${SORTED_UNIQUE_GENES_2K[COUNTER_ARR_2K]}" >> $FILE_2K
                                    let COUNTER_ARR_2K++
                                done
                                echo -e "" >> $FILE_2K    # insert a new line at the end of the motif data
                                
                                # counter for indices in arrays (SORTED_GENES_1K)
                                COUNTER_ARR_1K_ALL=0
                                # fill the *1k_all.txt file with data
                                echo -e "${GENE_SET_NAME_1K}" >> $FILE_1K_ALL   # insert the gene set name to the file (1k_all)
                                # append the genes to the assciated gene set name (1k_all)
                                for elem in "${SORTED_GENES_1K[@]}"; do
                                    echo -e "${SORTED_GENES_1K[COUNTER_ARR_1K_ALL]}" >> $FILE_1K_ALL
                                    let COUNTER_ARR_1K_ALL++
                                done
                                echo -e "" >> $FILE_1K_ALL    # insert a new line at the end of the motif data

                                # counter for indices in arrays (SORTED_GENES_2K)
                                COUNTER_ARR_2K_ALL=0
                                # fill the *2k_all.txt file with data
                                echo -e "${GENE_SET_NAME_2K}" >> $FILE_2K_ALL   # insert the gene set name to the file (2k_all)
                                # append the genes to the assciated gene set name (2k_all)
                                for elem in "${SORTED_GENES_2K[@]}"; do
                                    echo -e "${SORTED_GENES_1K[COUNTER_ARR_2K_ALL]}" >> $FILE_2K_ALL
                                    let COUNTER_ARR_2K_ALL++
                                done
                                echo -e "" >> $FILE_2K_ALL    # insert a new line at the end of the motif data
                            fi
                        else
                            echo "The motif $MOTIF of the cell type $CELL_TYPE of the tissue $TISSUE was already analyzed."
                        fi
                    done
                else
                    echo "The cell type $CELL_TYPE of the tissue $TISSUE was already analyzed."
                fi
            fi
        done
        echo "All cell types of the tissue $TISSUE have been analyzed."
    fi
done
echo "All tissues have been analyzed."
