#!/bin/bash

## script to check if the expected log files are existing and without errors
## quick finder for errors with the motif_discovery_pipeline

# get script path
SPATH=$(dirname $0)
# read in config
CONF="${SPATH}/../tglobal_vars.cnf"
while read LINE; do declare "$LINE"; done < $CONF

# input parameters
DIR="${PROJECT_DIR}/runs"
FILE_PATH=$DIR  # put into the evaluation folder

# output file
FILE_NAME=check_logs.txt
FILE=${FILE_PATH}/$FILE_NAME
if ! [ -f $FILE ]; then # check if the file is already existing
    cd $FILE_PATH
    touch $FILE_NAME
    echo "$FILE_NAME was created."
    echo -e "tissue\tcell_type\tlog_files\terrors\t" >> $FILE   # write headers in the empty file
fi

# check if any log files are existing and missing
for TISSUE in $DIR/*; do
    if [ -d $TISSUE ]; then
        TISSUE=$(echo $TISSUE | rev | cut -d'/' -f-1 | rev) # extract the tissue name
        TISSUE_AVAILABLE=$(cat $FILE | grep -c "${TISSUE}") # check if the tissue is already analyzed
        if [ "$TISSUE_AVAILABLE" = 0  ]; then
            # analyze the log files
            for CELL_TYPE in $DIR/$TISSUE/*; do
                if [ -d $CELL_TYPE/motif_discovery_pipeline ]; then
                    CELL_TYPE=$(echo $CELL_TYPE | rev | cut -d'/' -f-1 | rev)   # extract the cell type name
                    # initialize the counter for the log files and the default boolean value for "ERROR" as a String
                    COUNTER_LOGS=0
                    ERROR="no"
                    if [ -d $DIR/$TISSUE/$CELL_TYPE/motif_discovery_pipeline/logs ]; then
                        # count log files and check errors in files
                        for LOG_FILE in $DIR/$TISSUE/$CELL_TYPE/motif_discovery_pipeline/logs/*.log; do
                            let COUNTER_LOGS++
                            CHECK_CONTENT=$(cat $LOG_FILE | wc -l)  # check if the log files are empty
                            CHECK_ERROR=$(cat $LOG_FILE | grep -ic "Error") # check if the log files includes error messages
                            CHECK_EXCEPTION=$(cat $LOG_FILE | grep -ic "Exception") # check if the log files includes exception messages
                            let ERROR_CHECKER=CHECK_ERROR+CHECK_EXCEPTION
                            # check if any log file is empty
                            if [ $CHECK_CONTENT -eq 0 ]; then
                                ERROR="yes"
                            else
                                # check if there is an error or exception message
                                if [ $ERROR_CHECKER -gt 0 ]; then
                                    ERROR="yes"
                                fi
                            fi
                        done
                        # compare the counted to the expected log files
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
