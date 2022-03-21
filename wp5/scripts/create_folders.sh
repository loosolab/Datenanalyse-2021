#!/bin/bash

# script to create the needed folder structure

# input parameters
TISSUE=$1
CELL_TYPE=$2

# create tissue folder
if [ -d "/mnt/workspace_stud/allstud/wp5/runs/$TISSUE" ] ; then
    echo "The directory $TISSUE already exists. The directory wasn't created."
else
    mkdir /mnt/workspace_stud/allstud/wp5/runs/$TISSUE
fi

# create cell type folder
if [ -d "/mnt/workspace_stud/allstud/wp5/runs/$TISSUE/$CELL_TYPE" ] ; then
    echo "The directory $CELL_TYPE already exists. The directory wasn't created."
    echo "All folders already exist. The creation of the folders is canceled."
    exit 0
else
    mkdir /mnt/workspace_stud/allstud/wp5/runs/$TISSUE/$CELL_TYPE
fi

# create subfolders: motif_discovery_pipeline, annotation, differential_binding, evaluations
mkdir /mnt/workspace_stud/allstud/wp5/runs/$TISSUE/$CELL_TYPE/motif_discovery_pipeline
mkdir /mnt/workspace_stud/allstud/wp5/runs/$TISSUE/$CELL_TYPE/annotation
mkdir /mnt/workspace_stud/allstud/wp5/runs/$TISSUE/$CELL_TYPE/similarity

# check if the new folders were created successfully
if [ -d "/mnt/workspace_stud/allstud/wp5/runs/$TISSUE/$CELL_TYPE/evaluations" ] ; then
    echo "Folder structure for output files were created successfully."
else
    echo "Building folder structure failed."
    exit 1
fi
