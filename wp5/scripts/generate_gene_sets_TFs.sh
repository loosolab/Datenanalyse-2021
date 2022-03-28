#!/bin/bash

## script to generate txt files for the transcription factor (TF) gene set analyzation

# get script path
SPATH=$(dirname $0)
# read in config
CONF="${SPATH}/../global_vars.cnf"
while read LINE; do declare "$LINE"; done < $CONF

# input parameters
DIR="${TBSDIR}"
FILE_PATH=$DIR  # same as directory
FILE_PATH_SAVE="${PROJECT_DIR}/runs"   # path to save the files

# generate and fill txt files for the TF gene set analyzation
for TISSUE in $DIR/*; do
    if [ -d $TISSUE ]; then
        TISSUE=$(echo $TISSUE | rev | cut -d'/' -f-1 | rev) # extract the tissue name
        FILE_PATH_TISSUE="${FILE_PATH}/$TISSUE"
        FILE_PATH_SAVE_FILES="${FILE_PATH_SAVE}/$TISSUE"
        DT_TISSUE=$(date '+%d/%m/%Y %H:%M:%S')   # start time of the tissue
        echo "Starting ${TISSUE} at ${DT_TISSUE}."
        
        # generate the file for the distance up to 1000 bp with uniq genes
        FILE_NAME_1K=${TISSUE}_TFs_1k.txt
        FILE_PATH_SAVE_1K="$FILE_PATH_SAVE_FILES"   # path to save the file
        FILE_1K="${FILE_PATH_SAVE_1K}/$FILE_NAME_1K"
        if ! [ -f $FILE_1K ]; then  # only generate file if it isn't already existing
            cd $FILE_PATH_SAVE_1K
            touch $FILE_NAME_1K
            echo "$FILE_NAME_1K was created."
        fi
        # generate the file for the distance up to 2000 bp with unique genes
        FILE_NAME_2K=${TISSUE}_TFs_2k.txt
        FILE_PATH_SAVE_2K="$FILE_PATH_SAVE_FILES"   # path to save the file
        FILE_2K="${FILE_PATH_SAVE_2K}/$FILE_NAME_2K"
        if ! [ -f $FILE_2K ]; then  # only generate file if it isn't already existing
            cd $FILE_PATH_SAVE_2K
            touch $FILE_NAME_2K
            echo "$FILE_NAME_2K was created."
        fi
        # generate the file for the distance up to 1000 bp with all genes
        FILE_NAME_1K_ALL=${TISSUE}_TFs_1k_all.txt
        FILE_PATH_SAVE_1K_ALL="$FILE_PATH_SAVE_FILES"   # path to save the file
        FILE_1K_ALL="${FILE_PATH_SAVE_1K_ALL}/$FILE_NAME_1K_ALL"
        if ! [ -f $FILE_1K_ALL ]; then  # only generate file if it isn't already existing
            cd $FILE_PATH_SAVE_1K_ALL
            touch $FILE_NAME_1K_ALL
            echo "$FILE_NAME_1K_ALL was created."
        fi
        # generate the file for the distance up to 2000 bp with all genes
        FILE_NAME_2K_ALL=${TISSUE}_TFs_2k_all.txt
        FILE_PATH_SAVE_2K_ALL="$FILE_PATH_SAVE_FILES"   # path to save the file
        FILE_2K_ALL="${FILE_PATH_SAVE_2K_ALL}/$FILE_NAME_2K_ALL"
        if ! [ -f $FILE_2K_ALL ]; then  # only generate file if it isn't already existing
            cd $FILE_PATH_SAVE_2K_ALL
            touch $FILE_NAME_2K_ALL
            echo "$FILE_NAME_2K_ALL was created."
        fi
           
        FILE_PATH_ANNOTATION="${FILE_PATH_TISSUE}/output/TFBS"
                
        # check if the tissue is already analyzed
        TF_COUNTER=0  # counter for the to be analyzed TFs
        PUNCTUATION_CHECKER=$(cat $FILE_1K | grep -c "#")   # counts the analyzed TFs
        # counts the to be analyzed TFs
        if [ -d $FILE_PATH_ANNOTATION ]; then
            for TF in $FILE_PATH_ANNOTATION/*; do
                let TF_COUNTER++
            done
        fi
        # compare the counts
        if [ $TF_COUNTER -gt $PUNCTUATION_CHECKER ]; then
            for TF in $FILE_PATH_ANNOTATION/*; do
                # check if the file is a directory -> TF and nothing else
                if [ -d $TF ]; then
                    TF=$(echo $TF | rev | cut -d'/' -f-1 | rev)    # extract the TF name
                    # check if the TF has already been analyzed
                    TF_CHECKER=$(cat $FILE_1K | grep -c $TF)   #  check if the TF is in the file
                    if [ $TF_CHECKER -eq 0 ]; then
                        DT_TF=$(date '+%d/%m/%Y %H:%M:%S')   # start time of the TF
                        echo "Starting ${TF} at ${DT_TF}."
                                    
                        GENE_SET_NAME_1K="#${TF}_1k" # column 1 -> 1k
                        GENE_SET_NAME_2K="#${TF}_2k" # column 1 -> 2k
                        GENES_1K=()                  # columns after 2 -> 1k
                        GENES_2K=()                  # columns after 2 -> 2k
                                    
                        cd $FILE_PATH_ANNOTATION/$TF
                                        
                        # read *overview.txt
                        # check if the feature = gene $11, check the distance $16, add the gene name $20 to array
                        OVERVIEW=($(find . | grep overview.txt | sed 's,./,,g'))
                        GENES_1K=($(cat $OVERVIEW | awk -F'\t' '{if ($11 == "gene" && $16 != "NA" && $16 <= 1000 && $20 != "NA") print $20;}'))
                        GENES_2K=($(cat $OVERVIEW | awk -F'\t' '{if ($11 == "gene" && $16 != "NA" && $16 <= 2000 && $20 != "NA") print $20}'))
                                    
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
                    else
                        echo "The transcription factor $TF of the tissue $TISSUE was already analyzed."
                    fi
                fi
            done
        else
            echo "All transcription factors of the tissue $TISSUE have been analyzed."
        fi
    fi
done
echo "All tissues have been analyzed."
