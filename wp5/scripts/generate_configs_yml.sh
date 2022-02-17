#!/bin/bash

# input parameters 1 -> file path output directory
TISSUE=$1
CELL_TYPE=$2

# create new folders and subfolders
./create_folders.sh "$TISSUE" "$CELL_TYPE"

# output directory
OUTPUT_DIRECTORY="/mnt/workspace_stud/allstud/wp5/runs/$TISSUE/$CELL_TYPE/motif_discovery_pipeline"

# create new config file
cp /mnt/workspace_stud/allstud/wp5/source_files/config.yml /mnt/workspace_stud/allstud/wp5/configs/config_${TISSUE}_${CELL_TYPE}.yml

# used file for manipulation
FILE="/mnt/workspace_stud/allstud/wp5/configs/config_${TISSUE}_${CELL_TYPE}.yml"

# file paths
GENOME_FASTA="/mnt/workspace_stud/allstud/homo_sapiens.104.mainChr.fa"
SCORE_BIGWIG="/mnt/workspace_stud/allstud/wp3/tissues/$TISSUE/snakemakeout/footprinting/${CELL_TYPE}_footprints.bw"
PEAK_BED="/mnt/workspace_stud/allstud/wp3/tissues/$TISSUE/snakemakeout/peak_calling/ ${CELL_TYPE}_union.bed"
MOTIF_FILE="/mnt/workspace_stud/allstud/wp3/tissues/$TISSUE/snakemakeout/all_motifs.txt"

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
