#!/bin/bash

## script to create the needed folders for wp5 analysis

# get script path
SPATH=$(dirname $0)
# read in config
CONF="${SPATH}/./global_vars.cnf"
while read LINE; do declare "$LINE"; done < $CONF

# check if the folders are alredy existing
FILE_PATH="${PROJECT_DIR}"
FILE_PATH_CONFIGS="${PROJECT_DIR}/configs"
FILE_PATH_RUNS="${PROJECT_DIR}/runs"

# only generate file if the directories aren't already existing
if ! [ -d $FILE_PATH_CONFIGS ]; then
    if ! [ -d $FILE_PATH_RUNS ]; then
        # two needed folders for automated analyzation
        mkdir $PROJECT_DIR/configs       # store the produced config files
        mkdir $PROJECT_DIR/runs          # store the motif discovery output and another analyzation output
    fi
fi
