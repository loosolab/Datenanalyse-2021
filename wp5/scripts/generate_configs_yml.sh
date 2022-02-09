#!/bin/bash

# input parameters 1 -> file path output directory
TISSUE=$1
CELL_TYPE=$2

# TODO change output directory path
# output directory
OUTPUT_DIRECTORY="/mnt/workspace_stud/stud12/outputs_motif_discovery_pipeline/$TISSUE/$CELL_TYPE"

# create new config file
cp config.yml config_${TISSUE}_${CELL_TYPE}.yml

# TODO change path
# used file for manipulation
FILE="/mnt/workspace_stud/stud12/bash_scripts/config_${TISSUE}_${CELL_TYPE}.yml"

# TODO change paths
# input parameters 2 -> file paths
GENOME_FASTA="/mnt/workspace_stud/stud12/data_wp3/$TISSUE/flatfiles/$3"        # TODO set unified file and always call this?
SCORE_BIGWIG="/mnt/workspace_stud/stud12/data_wp3/$TISSUE/footprinting/$4"
PEAK_BED="/mnt/workspace_stud/stud12/data_wp3/$TISSUE/peak_calling/$5"
MOTIF_FILE="/mnt/workspace_stud/stud12/data_wp3/$TISSUE/motifs/$6"

# file manipulations
sed -i 's,^.*output:.*$,'"  output: \'$OUTPUT_DIRECTORY\'"',' $FILE
sed -i 's,^.*genome_fasta:.*$,'"  genome_fasta: \'$GENOME_FASTA\'"',' $FILE
sed -i 's,^.*score_bigwig:.*$,'"  score_bigwig: \'$SCORE_BIGWIG\'"',' $FILE
sed -i 's,^.*peak_bed:.*$,'"  peak_bed: \'$PEAK_BED\'"',' $FILE
sed -i 's,^.*motif_file:.*$,'"  motif_file: \'$MOTIF_FILE\'"',' $FILE

# check if new output directory is inserted into file
OUTPUT_CHECKER=$(grep -c "$OUTPUT_DIRECTORY" $FILE)

# check if the new config file was created successfully
if [ -f "config_${TISSUE}_${CELL_TYPE}.yml" ] ; then
    if [ $OUTPUT_CHECKER = 1 ] ; then
        echo "config_${TISSUE}_${CELL_TYPE}.yml was created successfully"
    else
        echo "config_${TISSUE}_${CELL_TYPE}.yml wasn't created successfully"
    fi
fi
