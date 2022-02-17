#!/bin/bash

# input parameters
TISSUE=$1
CELL_TYPE=$2

# create new folders and subfolders
./create_folders.sh "$TISSUE" "$CELL_TYPE"

# create new config file
cp /mnt/workspace_stud/allstud/wp5/configs/config_${TISSUE}_${CELL_TYPE}.yml /mnt/workspace_stud/allstud/wp5/configs/config_${TISSUE}_${CELL_TYPE}_with_annotation.yml

# choose gtf file
GTF="/mnt/workspace_stud/allstud/homo_sapiens.104.mainChr.gtf"

# choose config file
FILE="/mnt/workspace_stud/allstud/wp5/configs/config_${TISSUE}_${CELL_TYPE}_with_annotation.yml"

# choose UROPA file
UROPA="/mnt/workspace_stud/allstud/wp5/source_files/uropa_template_${TISSUE}_${CELL_TYPE}.json"

# output directory
OUTPUT_DIRECTORY="/mnt/workspace_stud/allstud/wp5/runs/$TISSUE/$CELL_TYPE/annotation"

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
