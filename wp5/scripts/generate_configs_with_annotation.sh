#!/bin/bash

# script to generate config files with an annotation step

# get script path
SPATH=$(dirname $0)
# read in config
CONF="${SPATH}/../test.conf"
while read LINE; do declare "$LINE"; done < $CONF

# input parameters
TISSUE=$1
CELL_TYPE=$2

# create new folders and subfolders
./create_folders.sh "$TISSUE" "$CELL_TYPE"

# create new config file
if [ -f "${PROJECT_DIR}/configs/config_${TISSUE}_${CELL_TYPE}_with_annotation.yml" ] ; then
    echo "The file config_${TISSUE}_${CELL_TYPE}_with_annotation.yml already exists. The file wasn't created."
    exit 0
else
    cp ${PROJECT_DIR}/configs/config_${TISSUE}_${CELL_TYPE}.yml ${PROJECT_DIR}/configs/config_${TISSUE}_${CELL_TYPE}_with_annotation.yml
fi

# choose gtf file TODO: change storing place of gtf
GTF="${PROJECT_DIR}/../homo_sapiens.104.mainChr.gtf"

# choose config file
FILE="${PROJECT_DIR}/configs/config_${TISSUE}_${CELL_TYPE}_with_annotation.yml"

# choose UROPA file
UROPA="${PROJECT_DIR}/source_files/uropa_template_${TISSUE}_${CELL_TYPE}.json"

# output directory
OUTPUT_DIRECTORY="${PROJECT_DIR}/runs/$TISSUE/$CELL_TYPE/annotation"

# file manipulation
sed -i 's,^.*output:.*$,'"  output: \'$OUTPUT_DIRECTORY\'"',' $FILE
sed -i 's,^.*gtf:.*$, '"  gtf: \'$GTF\'"',' $FILE
sed -i 's,^.*config_template:.*$, '"  config_template: \'$UROPA\'"',' $FILE

# check if new output directory is inserted into file
OUTPUT_CHECKER=$(grep -c "$OUTPUT_DIRECTORY" $FILE)

# check if the new config file was created successfully
if [ -f "config_${TISSUE}_${CELL_TYPE}.yml" ] ; then
    if [ $OUTPUT_CHECKER = 1 ] ; then
        echo "config_${TISSUE}_${CELL_TYPE}.yml was created successfully"
    else
        echo "Creating config file failed."
        exit 1
    fi
fi
