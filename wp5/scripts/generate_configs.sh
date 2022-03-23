#!/bin/bash
## Script to generate config for each celltype of each tissue

# Get path to this script
SPATH=$(dirname $0)
# read in config
CONF="${SPATH}/../test.conf"
while read LINE; do declare "$LINE"; done < $CONF

# generate config for each celltype of each tissue
for TISSUE in $TBSDIR/*/; do
    TIS_NAME=$(basename $TISSUE)
    for SUB_TYPE in ${TISSUE}output/footprinting/*_footprints.bw; do
	    SUB_NAME=$(basename "$SUB_TYPE" "_footprints.bw")
	    if [[ ${ANN_CHECKER} == "yes" ]]; then
            ${SPATH}/utils/generate_configs_with_annotation.sh $TIS_NAME $SUB_NAME
	    else
            ${SPATH}/utils/generate_configs_yml.sh $TIS_NAME $SUB_NAME
	    fi
    done    
done
