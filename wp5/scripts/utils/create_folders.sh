#!/bin/bash

## script to create the needed folder structure

# get script path
SPATH=$(dirname $0)
# read in config
CONF="${SPATH}/../../tglobal_vars.cnf"
while read LINE; do declare "$LINE"; done < $CONF

# input parameters
TISSUE=$1
CELL_TYPE=$2

# create tissue folder
if [ -d "${PROJECT_DIR}/runs/$TISSUE" ] ; then
    echo "The directory $TISSUE already exists. The directory wasn't created."
else
    mkdir ${PROJECT_DIR}/runs/$TISSUE
fi

# create cell type folder
if [ -d "${PROJECT_DIR}/runs/$TISSUE/$CELL_TYPE" ] ; then
    echo "The directory $CELL_TYPE already exists. The directory wasn't created."
    echo "All folders already exist. The creation of the folders is canceled."
    exit 0
else
    mkdir ${PROJECT_DIR}/runs/$TISSUE/$CELL_TYPE
fi

# create subfolders: motif_discovery_pipeline, annotation, differential_binding, evaluations
mkdir ${PROJECT_DIR}/runs/$TISSUE/$CELL_TYPE/motif_discovery_pipeline
mkdir ${PROJECT_DIR}/runs/$TISSUE/$CELL_TYPE/annotation
mkdir ${PROJECT_DIR}/runs/$TISSUE/$CELL_TYPE/similarity

# check if the new folders were created successfully
if [ -d "${PROJECT_DIR}/runs/$TISSUE/$CELL_TYPE/annotation" ] ; then
    echo "Folder structure for output files were created successfully."
else
    echo "Building folder structure failed."
    exit 1
fi
