#!/bin/bash
## read input (Tissues and cells of the two conditions + output directory) 
T1=$1
T2=$2
C1=$3
C2=$4
OUTPUT_DIRECTORY="${5}/TOBIAS"

## define needed paths
P1="${T1}/${C1}/motif_discovery_pipeline/3_evaluation/motifs.meme"
P2="/mnt/workspace_stud/allstud/wp5/runs/${T2}/${C2}/motif_discovery_pipeline/3_evaluation/motifs.meme"
SIGNALS1="../test_data_tobias_output/footprinting/mDuxPos_footprints.bw"
SIGNALS2="../test_data_tobias_output/footprinting/mDuxNeg_footprints.bw"
GENOME="../test_data_tobias_output/flatfiles/mus_musculus.94.mainChr.fa"
PEAKS="../test_data_tobias_output/peak_calling/all_merged.bed"

## differential binding w tobias:
# activate TOBIAS ENV
source /opt/miniconda/bin/activate TOBIAS_ENV
# run
TOBIAS BINDetect --motifs $P1 $P2  --signals $SIGNALS1 $SIGNALS2   --genome $GENOME --peaks $PEAKS --cond_names "${T1}_${C1}" "${T2}_${C2}"  --outdir "${OUTPUT_DIRECTORY}"
