#!/bin/bash

## Script to generate config for each celltype of each tissue

# get script path
SPATH=$(dirname $0)
# read in config
CONF="${SPATH}/../global_vars.cnf"
while read LINE; do declare "$LINE"; done < $CONF

# generate config for each celltype of each tissue
for TISSUE in $TBSDIR/*/; do
    TIS_NAME=$(basename $TISSUE)
    for SUB_TYPE in ${TISSUE}output/footprinting/*_footprints.bw; do
	    SUB_NAME=$(basename "$SUB_TYPE" "_footprints.bw")
	    if [ $ANN_CHECKER = "yes" ]; then
            # create new folders and subfolders
            ${SPATH}/utils/create_folders.sh $TIS_NAME $SUB_NAME
            ${SPATH}/utils/generate_configs_with_annotation.sh $TIS_NAME $SUB_NAME $GENOME $GTF
	    else
            # create new folders and subfolders
            ${SPATH}/utils/create_folders.sh $TIS_NAME $SUB_NAME
            ${SPATH}/utils/generate_configs_yml.sh $TIS_NAME $SUB_NAME $GENOME 
	    fi
    done    
done
