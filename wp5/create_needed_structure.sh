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
FILE_PATH_SOURCE_FILES="${PROJECT_DIR}/source_files"

# only generate file if the directories aren't already existing
if ! [ -d $FILE_PATH_CONFIGS ]; then
    if ! [ -d $FILE_PATH_RUNS ]; then
        if ! [ -d $FILE_PATH_SOURCE_FILES ]; then
            # three needed folders for automated analyzation
            mkdir $PROJECT_DIR/configs       # store the produced config files
            mkdir $PROJECT_DIR/runs          # store the motif discovery output and another analyzation output
            mkdir $PROJECT_DIR/source_files  # store default config.yml, uropa_template.json file, genome.fa, genome.gtf
        fi
    fi
fi

# generate info txt -> how to prepare the source_files folder to start the analysis
FILE_NAME="prepare_source_files.txt"
FILE_PATH_SAVE="$FILE_PATH_SOURCE_FILES"   # path to save the file
FILE="${FILE_PATH_SAVE}/$FILE_NAME"
if ! [ -f $FILE ]; then  # only generate file if it isn't already existing
    cd $FILE_PATH_SAVE
    touch $FILE_NAME
    echo "Please copy the following four files to that directory:" >> $FILE
    echo "1: Please copy the default config.yml from the motif-discovery-pipeline." >> $FILE
    echo "2: Please copy the default uropa_template.json from the motif-discovery-pipeline." >> $FILE
    echo "3: Please copy the genome.fa you want to use for your analysis. If you use the output from wp3 please copy the genome.fa from one of the outputs (in the subfolder wp3/new_tissiues/some_tissue/output/flatfiles) of wp3." >> $FILE
    echo "4: Please copy the genome.gtf you want to use for your analysis. If you use the output from wp3 please copy the genome.gtf from one of the outputs (in the subfolder wp3/new_tissiues/some_tissue/output/flatfiles) of wp3." >> $FILE
    echo "$FILE_NAME was created."
fi
