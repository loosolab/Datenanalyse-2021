#!/bin/bash

T1=$1
T2=$2
C1=$3
C2=$4
OUT_T="${5}/TOBIAS_${T1}_${C1}_${T2}_${C2}/"

P1="/mnt/workspace_stud/allstud/wp5/runs/${T1}/${C1}/motif_discovery_pipeline/3_evaluation/motifs.meme"
P2="/mnt/workspace_stud/allstud/wp5/runs/${T2}/${C2}/motif_discovery_pipeline/3_evaluation/motifs.meme"
MOTIF_DB="" # TODO ?
SIGNALS1="../test_data_tobias_output/footprinting/mDuxPos_footprints.bw"
SIGNALS2="../test_data_tobias_output/footprinting/mDuxNeg_footprints.bw"
GENOME="../test_data_tobias_output/flatfiles/mus_musculus.94.mainChr.fa"
PEAKS="../test_data_tobias_output/peak_calling/all_merged.bed"
## differential binding w tobias:
TOBIAS BINDetect --motifs $P1 $P2 $MOTIF_DB  --signals $SIGNALS1 $SIGNALS2   --genome $GENOME --peaks $PEAKS --cond_names "${T1}_${C1}" "${T2}_${C2}"  --outdir $OUT_T

# TODO
# join new motifs with db 

# make new configs
    # get config
    # extract names
    # copy config
    # enter names
    # change paths

# start runs 

# extract tsvs

# evaluate #TODO edit params
/mnt/workspace_stud/allstud/wp5/scripts/evalDiff.py CONDITION1.tsv CONDITION2.tsv --BD-comp "${OUT_T}BINDetect_results.txt" --out "${T1}_${C1}_${T2}_${C2}_Eval" --get-outliers --cond_names "${T1}_${C1}" "${T2}_${C2}" 