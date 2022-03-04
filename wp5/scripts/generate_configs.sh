#!/bin/bash
# TOBIAS output directory
# Path to tobias directory
TBSDIR="/mnt/workspace_stud/allstud/wp3/tissues"

for TISSUE in $TBSDIR/*/; do
    TIS_NAME=$(basename $TISSUE)
    for SUB_TYPE in ${TISSUE}snakemakeout/footprinting/*_footprints.bw; do
	    SUB_NAME=$(basename "$SUB_TYPE" "_footprints.bw")
	    ./generate_configs_yml.sh $TIS_NAME $SUB_NAME
    done    
done