#!/bin/bash

## script to create the needed folders for wp5 analysis

# get script path
SPATH=$(dirname $0)
# read in config
CONF="${SPATH}/../../global_vars.cnf"
while read LINE; do declare "$LINE"; done < $CONF

# paths of needed folders for automated analyzation
FILE_PATH_CONFIGS="${PROJECT_DIR}/configs"
FILE_PATH_RUNS="${PROJECT_DIR}/runs"
FILE_PATH_SIM="${PROJECT_DIR}/similarity"

# only generate directories if they do not already exist
if ! [ -d $FILE_PATH_CONFIGS ]; then
    mkdir $FILE_PATH_CONFIGS       # store the produced config files
fi
if ! [ -d $FILE_PATH_RUNS ]; then
    mkdir $FILE_PATH_RUNS         # store the motif discovery output and another analyzation output
fi
if ! [ -d $FILE_PATH_SIM ]; then
    mkdir $FILE_PATH_SIM         # store the motiv similarity evaluation output
fi