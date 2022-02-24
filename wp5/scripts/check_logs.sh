#!/bin/bash

# input parameters
DIR="/mnt/workspace_stud/allstud/wp5/runs"
# LOG1="footprint_extraction.log"
# LOG2="meme.log"
# LOG3="motif_discovery.log"
# LOG4="rescan.log"
FILE_PATH="/mnt/workspace_stud/allstud/wp5/runs"  # put into the evaluation folder

# output file
FILE_NAME=check_logs.txt
FILE=${FILE_PATH}/$FILE_NAME
if ! [ -f $FILE ]; then
    cd $FILE_PATH
    touch $FILE_NAME
    echo "$FILE_NAME was created."
    echo -e "tissue\tcell_type\tlog_files\terrors\t" >> $FILE
fi

# check if any log files are existing and missing
for TISSUE in $DIR/*; do
    if [ -d $TISSUE ]; then
        TISSUE=$(echo $TISSUE | rev | cut -d'/' -f-1 | rev)
        TISSUE_AVAILABLE=$(cat $FILE | grep -c "${TISSUE}")
        if [ "$TISSUE_AVAILABLE" = 0  ]; then
            for CELL_TYPE in $DIR/$TISSUE/*; do
                if [ -d $CELL_TYPE/motif_discovery_pipeline ]; then
                    CELL_TYPE=$(echo $CELL_TYPE | rev | cut -d'/' -f-1 | rev)
                    COUNTER_LOGS=0
                    ERROR="no"
                    if [ -d $DIR/$TISSUE/$CELL_TYPE/motif_discovery_pipeline/logs ]; then
                        for LOG_FILE in $DIR/$TISSUE/$CELL_TYPE/motif_discovery_pipeline/logs/*.log; do
                            let COUNTER_LOGS++
                            CHECK_CONTENT=$(cat $LOG_FILE | wc -l)
                            CHECK_ERROR=$(cat $LOG_FILE | grep -ic "Error")
                            CHECK_EXCEPTION=$(cat $LOG_FILE | grep -ic "Exception")
                            let ERROR_CHECKER=CHECK_ERROR+CHECK_EXCEPTION
                            if [ $CHECK_CONTENT -eq 0 ]; then
                                ERROR="yes"
                            else
                                if [ $ERROR_CHECKER -gt 0 ]; then
                                    ERROR="yes"
                                fi
                            fi
                        done
                        if [ $COUNTER_LOGS -lt 4 ]; then
                            ERROR="yes"
                        fi
                    fi
                    echo -e "${TISSUE}\t${CELL_TYPE}\t${COUNTER_LOGS}\t${ERROR}\t" >> $FILE
                fi
            done
        else
            echo "$TISSUE is already analyzed."
        fi
    fi
done
