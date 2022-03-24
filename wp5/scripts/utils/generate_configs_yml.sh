#!/bin/bash

## script to generate config files without an annotation step

# get script path
SPATH=$(dirname $0)
# read in config
CONF="${SPATH}/../../global_vars.cnf"
while read LINE; do declare "$LINE"; done < $CONF

# input parameters 1 -> file path output directory
TISSUE=$1
CELL_TYPE=$2

# output directory
OUTPUT_DIRECTORY="${PROJECT_DIR}/runs/$TISSUE/$CELL_TYPE/motif_discovery_pipeline"

# create new config file
if [ -f "${PROJECT_DIR}/configs/config_${TISSUE}_${CELL_TYPE}.yml" ] ; then
    echo "The file config_${TISSUE}_${CELL_TYPE}.yml already exists. The file wasn't created."
    exit 0
else
    cp $PROJECT_DIR/source_files/config.yml ${PROJECT_DIR}/configs/config_${TISSUE}_${CELL_TYPE}.yml
fi

# used file for manipulation
FILE="${PROJECT_DIR}/configs/config_${TISSUE}_${CELL_TYPE}.yml"

# file paths TODO: change storing place of genome fasta
GENOME_FASTA="${PROJECT_DIR}/../homo_sapiens.104.mainChr.fa"
SCORE_BIGWIG="${TBSDIR}/$TISSUE/output/footprinting/${CELL_TYPE}_footprints.bw"
PEAK_BED="${TBSDIR}/$TISSUE/output/peak_calling/${CELL_TYPE}_union.bed"
MOTIF_FILE="${TBSDIR}/$TISSUE/output/motifs/all_motifs.txt"

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
        echo "Creating config file failed."
        exit 1
    fi
fi
