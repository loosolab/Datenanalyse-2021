#!/bin/bash
## Script to generate config for each celltype of each tissue

# TOBIAS output directory
# Path to tobias directory ; TODO
TBSDIR="/mnt/workspace_stud/allstud/wp3/new_tissiues"
# Get path to this script
SPATH=$(dirname $0)

# generate config for each celltype of each tissue
for TISSUE in $TBSDIR/*/; do
    TIS_NAME=$(basename $TISSUE)
    for SUB_TYPE in ${TISSUE}snakemakeout/footprinting/*_footprints.bw; do
	    SUB_NAME=$(basename "$SUB_TYPE" "_footprints.bw")
	    $(SPATH)/utils/generate_configs_yml.sh $TIS_NAME $SUB_NAME
    done    
done
