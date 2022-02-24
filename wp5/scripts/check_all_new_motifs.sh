#!/bin/bash

# input parameters
DIR="/mnt/workspace_stud/allstud/wp5/runs"
FILE_PATH="/mnt/workspace_stud/allstud/wp5/runs"

# output: all_motifs file
FILE_NAME=check_new_motifs_all.txt
FILE=${FILE_PATH}/$FILE_NAME
if ! [ -f $FILE ]; then
    cd $FILE_PATH
    touch $FILE_NAME
    echo "$FILE_NAME was created."
    echo -e "tissue\tcell_type\tsequence\tmotif_name\t" >> $FILE 

    # check all check_new_motifs files and concatenate them
    for TISSUE in $DIR/*; do
        if [ -d $TISSUE ]; then
            TISSUE=$(echo $TISSUE | rev | cut -d'/' -f-1 | rev)
            # check for all check_new_motif files
            FILE_PATH_TISSUE="${FILE_PATH}/$TISSUE"
            FILE_NAME_TISSUE=check_new_motifs_${TISSUE}.txt
            
            cd $FILE_PATH_TISSUE
            # variable to ignore the first input line
            IGNORER=0

            # read check_new_motifs_tissue.txt files
            while read line
            do
                if [ $IGNORER -eq 1 ]; then
                # strip the input line and write the values into an array
                stripped_line=($(echo $line | tr "\t" "\n"))
                # assign the values to variables
                CURRENT_TISSUE=$(echo ${stripped_line[0]})
                CURRENT_CELL_TYPE=$(echo ${stripped_line[1]})
                ID=$(echo ${stripped_line[2]})
                CONSENUS=$(echo ${stripped_line[3]})                
                # add line to file
                    echo -e "${CURRENT_TISSUE}\t${CURRENT_CELL_TYPE}\t${ID}\t${CONSENUS}\t" >> $FILE
                else
                # start loop after the first input line
                    let IGNORER++
                fi
            done < $FILE_NAME_TISSUE
        fi
    done

else
    echo "$FILE_NAME was already created."
fi
