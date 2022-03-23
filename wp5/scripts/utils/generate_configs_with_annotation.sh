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

if [[ ${ANN_CHECKER} == "yes" ]]; then
    # create new config file
    if [ -f "${PROJECT_DIR}/configs/config_${TISSUE}_${CELL_TYPE}_with_annotation.yml" ] ; then
        echo "The file config_${TISSUE}_${CELL_TYPE}_with_annotation.yml already exists. The file wasn't created."
        exit 0
    else
        cp ${PROJECT_DIR}/source_files/config.yml ${PROJECT_DIR}/configs/config_${TISSUE}_${CELL_TYPE}_with_annotation.yml
    fi

    # used file for manipulation
    FILE="${PROJECT_DIR}/configs/config_${TISSUE}_${CELL_TYPE}_with_annotation.yml"

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
else
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
fi
