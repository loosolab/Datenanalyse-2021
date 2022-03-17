#!/bin/bash

## read input (Tissues and cells of the two conditions + output directory) 
T1=$1
T2=$2
C1=$3
C2=$4
OUTPUT_DIRECTORY="${5}/MD_Pipeline"

## define needed paths
MOTIF_FILE="${OUTPUT_DIRECTORY}/joined_motifs.jaspar"
MOTIFS_C1="${T1}/${C1}/motif_discovery_pipeline/3_evaluation/motifs.meme"
MOTIFS_C2="/mnt/workspace_stud/allstud/wp5/runs/${T2}/${C2}/motif_discovery_pipeline/3_evaluation/motifs.meme"

## combine motifs:
# activate TOBIAS ENV
source /opt/miniconda/bin/activate TOBIAS_ENV
# join motifs
TOBIAS FormatMotifs --input $MOTIFS_C1 $MOTIFS_C2 --task join --output $MOTIF_FILE
# deactivate env
conda deactivate


## make new configs
# copy config
CONF1_PATH="${OUTPUT_DIRECTORY}/configs/${T2}_${C2}_under_${T1}_${C1}.yml"
CONF2_PATH="${OUTPUT_DIRECTORY}/configs/${T1}_${C1}_under_${T2}_${C2}.yml"
cp "/mnt/workspace_stud/allstud/wp5/configs/config_${T1}_${C1}.yml" $CONF1_PATH
cp "/mnt/workspace_stud/allstud/wp5/configs/config_${T2}_${C2}.yml" $CONF2_PATH

# file manipulations
sed -i 's,^.*output:.*$,'"  output: \'${OUTPUT_DIRECTORY}/${T2}_${C2}_under_${T1}_${C1}/\'"',' $CONF1_PATH
sed -i 's,^.*motif_file:.*$,'"  motif_file: \'$MOTIF_FILE\'"',' $CONF1_PATH

sed -i 's,^.*output:.*$,'"  output: \'${OUTPUT_DIRECTORY}/${T1}_${C1}_under_${T2}_${C2}/\'"',' $CONF2_PATH
sed -i 's,^.*motif_file:.*$,'"  motif_file: \'$MOTIF_FILE\'"',' $CONF2_PATH

# extract and enter motif names in yaml
grep ">" $MOTIF_FILE | cut -f1 | while read motif_name; do sed -i 's,^    motifs:.*$,'"    motifs:\n      - \"$motif_name\""',' $CONF1_PATH;done
grep ">" $MOTIF_FILE | cut -f1 | while read motif_name; do sed -i 's,^    motifs:.*$,'"    motifs:\n      - \"$motif_name\""',' $CONF2_PATH;done

## start runs 
# activate conda
source /opt/miniconda/bin/activate MEME_ENV
# run motif discovery pipeline
./run_pipeline.sh "${OUTPUT_DIRECTORY}/configs"
# deactivate conda
conda deactivate