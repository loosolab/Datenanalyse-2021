#!/bin/bash
# TOBIAS output directory
# Path to tobias directory
TBSDIR = $1
# output File to write config files to 
CONFDIR = $2
# Directory where the pipeline results should be wrote to
OUTDIR = $3
# path to config template.
# TODO
CONFF = "../motif-discovery-pipeline/config.yml"

# path to genome and motif file
GENOME = ""
MOTIFDB = "${TBSDIR}/motifs/all_motifs.txt"
for filename in $TBSDIR/peak_calling/*_union.bed; do
    echo $filename
    # get pattern for naming
    NAME = "test"
    YMLNAME = $CONFDIR/config_$NAME.yml
    # copy yml
    #cp $CONFF $YMLNAME
    # set output directory
    sed '/output_dir/${OUTDIR}/' $YMLNAME
    # set genome file
    sed "/  genome_fasta: ''/  genome_fasta: '${GENOME}'/"
    # set score bigwig in /footprinting/ ..._footprinting.bw
    sed "/  score_bigwig: ''/  score_bigwig: '${TBSDIR}/footprinting/${NAME}_footprinting.bw'/"
    # set peak bed in /peak_calling/ ..._union.bed
    sed "/  peak_bed: ''/  peak_bed: '${TBSDIR}/peak_calling/${filename}'/"
    # set motif database
    sed "/  motif_file: ''/  motif_file: '${MOTIFDB}'/"

    # TODO start pipeline run for config?

done
