#!/bin/bash

# script to generate gmt files for the GSEA

# input parameters
DIR="/mnt/workspace_stud/allstud/wp5/runs"
FILE_PATH=$DIR  # same as directory

# generate and fill gmt files for the GSEA analyzation
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
                # generate the file for the distance up to 1000 bp
                FILE_NAME_1K=${TISSUE}_${CELL_TYPE}_1k.gmt
                FILE_PATH_SAVE_1K="${FILE_PATH_CELL_TYPE}/annotation"   # path to save the file
                FILE_1K="${FILE_PATH_SAVE_1K}/$FILE_NAME_1K"
                if ! [ -f $FILE_1K ]; then  # only generate file if it isn't already existing
                    cd $FILE_PATH_SAVE_1K
                    touch $FILE_NAME_1K
                    echo "$FILE_NAME_1K was created."
                fi
                # generate the file for the distance up to 2000 bp
                FILE_NAME_2K=${TISSUE}_${CELL_TYPE}_2k.gmt
                FILE_PATH_SAVE_2K="${FILE_PATH_CELL_TYPE}/annotation"   # path to save the file
                FILE_2K="${FILE_PATH_SAVE_2K}/$FILE_NAME_2K"
                if ! [ -f $FILE_2K ]; then  # only generate file if it isn't already existing
                    cd $FILE_PATH_SAVE_2K
                    touch $FILE_NAME_2K
                    echo "$FILE_NAME_2K was created."
                fi
                
                FILE_PATH_ANNOTATION="${FILE_PATH_CELL_TYPE}/motif_discovery_pipeline/4_annotation/uropa"
                
                # check if the cell types are already analyzed
                MOTIF_COUNTER=0  # counter for the to be analyzed motifs
                LINE_CHECKER=$(cat $FILE_1K | wc -l)   # counts the analyzed motifs
                # counts the to be analyzed motifs
                if [ -d $FILE_PATH_ANNOTATION ]; then
                    for MOTIF in $FILE_PATH_ANNOTATION/*; do
                        if [[ $MOTIF != "config" ]];then
                            let MOTIF_COUNTER++
                        fi
                    done
                fi
                # compare the counts
                if [ $MOTIF_COUNTER -ge $LINE_CHECKER ]; then
                    DT_CELL_TYPE=$(date '+%d/%m/%Y %H:%M:%S')   # start time of the cell type
                    echo "Starting ${CELL_TYPE} at ${DT_CELL_TYPE}."
                    for MOTIF in $FILE_PATH_ANNOTATION/*; do
                        MOTIF=$(echo $MOTIF | rev | cut -d'/' -f-1 | rev)    # extract the motif name
                        if [[ $MOTIF != "config" ]]; then
                            DT_MOTIF=$(date '+%d/%m/%Y %H:%M:%S')   # start time of the motif
                            echo "Starting ${MOTIF} at ${DT_MOTIF}."
                            
                            GENE_SET_NAME_1K="$MOTIF_1k"    # column 1 -> 1k
                            GENE_SET_NAME_2K="$MOTIF_2k"    # column 1 -> 2k
                            DESCRIPTION="na"                # column 2 -> 1k + 2k
                            GENES_1K=()                     # columns after 2 -> 1k
                            GENES_2K=()                     # columns after 2 -> 2k
                                
                            cd $FILE_PATH_ANNOTATION/$MOTIF
                            # variable to ignore the first input line
                            IGNORER=0
                                
                            # read *allhits.txt
                            while read line
                            do
                                if [ $IGNORER -eq 1 ]; then
                                    # strip the input line and write the values into an array
                                    stripped_line=($(echo $line | tr "\t" "\n"))
                                    # check the distance
                                    DISTANCE=(${stripped_line[11]})
                                    # check if the value for distance isn't "NA"
                                    if [[ $DISTANCE != "NA" ]]; then
                                        # check if the distance is less equal than 1k
                                        if [ $DISTANCE -le 1000 ]; then
                                            # assign the needed values to variables -> 1k
                                            GEN_1K=$(echo ${stripped_line[17]})
                                            # append the variable values to the arrays
                                            GENES_1K+=($GEN_1K)
                                        fi
                                        # check if the distance is less equal than 2k
                                        if [ $DISTANCE -le 2000 ]; then
                                            # assign the needed values to variables -> 2k
                                            GEN_2K=$(echo ${stripped_line[17]})
                                            # append the variable values to the arrays
                                            GENES_2K+=($GEN_2K)
                                        fi
                                    fi
                                else
                                    # start loop after the first input line
                                    let IGNORER++
                                fi
                            done < *allhits.txt
                                
                            # only choose the unique genes
                            SORTED_UNIQUE_GENES_1K=($(echo "${GENES_1K[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
                            SORTED_UNIQUE_GENES_2K=($(echo "${GENES_2K[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
                                                            
                            # counter for indices in arrays (GENES_1K)
                            COUNTER_ARR_1K=0
                            # fill the *1k.gmt file with data
                            echo -e "${GENE_SET_NAME_1K}" >> $FILE_1K   # insert the gene set name to the file (1k)
                            sed -i '/^'$GENE_SET_NAME_1K'/s!$!\t'$DESCRIPTION'!' $FILE_1K   # append the description to the assciated gene set name (1k)
                            # append the genes to the assciated gene set name (1k)
                            for elem in "${SORTED_UNIQUE_GENES_1K[@]}"; do
                                sed -i '/^'$GENE_SET_NAME_1K'/s!$!\t'${SORTED_UNIQUE_GENES_1K[COUNTER_ARR_1K]}'!' $FILE_1K
                                let COUNTER_ARR_1K++
                            done
                            sed -i '/^'$GENE_SET_NAME_1K'/s!$!\t!' $FILE_1K # insert a tab at the end of the line

                            # counter for indices in arrays (GENES_2K)
                            COUNTER_ARR_2K=0
                            # fill the *2k.gmt file with data
                            echo -e "${GENE_SET_NAME_2K}" >> $FILE_2K   # insert the gene set name to the file (2k)
                            sed -i '/^'$GENE_SET_NAME_2K'/s!$!\t'$DESCRIPTION'!' $FILE_2K   # append the description to the assciated gene set name (2k)
                            # append the genes to the assciated gene set name (2k)
                            for elem in "${SORTED_UNIQUE_GENES_2K[@]}"; do
                                sed -i '/^'$GENE_SET_NAME_2K'/s!$!\t'${SORTED_UNIQUE_GENES_2K[COUNTER_ARR_2K]}'!' $FILE_2K
                                let COUNTER_ARR_2K++
                            done
                            sed -i '/^'$GENE_SET_NAME_2K'/s!$!\t!' $FILE_2K # insert a tab at the end of the line
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
