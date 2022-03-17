#!/bin/bash

# remove last slash of directory path if existing
if [[ $1 == */ ]]; then
  OUT_DIR_PRE=${1::-1};
else
  OUT_DIR_PRE=$1; 
fi

echo "start"
DIR="/mnt/workspace_stud/allstud/wp5/runs/"
# get Tissue names
TISSUES=( $( find ${DIR} -mindepth 1 -maxdepth 1 -type d -printf '%f\n' ) )

# LOOP over tissues
START=0
END=${#TISSUES[@]}
# TODO: differential binding between tissues
#END2=${#TISSUES[@]}
#END=$(($END2 - 1))
 
for (( i=$START; i<$END; i++ ))
do
    #START2=$(($i + 1))
    T1=${TISSUES[$i]}
    CELL_TYPES=( $( find ${DIR}/${T1}/ -mindepth 1 -maxdepth 1 -type d -printf '%f\n' ) )
    END_C2=${#CELL_TYPES[@]}
    END_C1=$(($END_C2 - 1))#
    # iterate over all celltype combinations
    for (( j=$START; j<$END_C1; j++ )); do
        START_C2=$(($j + 1))
        for (( k=$START_C2; k<$END_C2; k++ )); do
            C1=${CELL_TYPES[$j]}
            C2=${CELL_TYPES[$k]} 
            OUT_DIR="${OUT_DIR_PRE}/DiffBind_${T1}_${C1}_${T1}_${C2}"
            # using TOBIAS
            ./DifferentialBinding_TOBIAS.sh $T1 $T1 ${C1} ${C2} $OUT_DIR
            # using Motif discovery pipeline
            ./DifferentialBinding_MDP.sh $T1 $T1 ${C1} ${C2} $OUT_DIR
            # extract tsvs
            CONDITION1_tsv="${OUTPUT_DIRECTORY}/MD_Pipeline/${T1}_${C2}_under_${T1}_${C1}/3_evaluation/motif_ranks.tsv"
            CONDITION2_tsv="${OUTPUT_DIRECTORY}/MD_Pipeline/${T1}_${C2}_under_${T1}_${C1}/3_evaluation/motif_ranks.tsv"
            # evaluate 
            /mnt/workspace_stud/allstud/wp5/scripts/evalDiff.py $CONDITION1_tsv $CONDITION2_tsv --BD-comp "${OUT_DIR}/TOBIAS/BINDetect_results.txt" --out "${OUT_DIR}/${T1}_${C1}_${T1}_${C2}_Eval" --get-outliers --cond_names "${T1}_${C1}" "${T1}_${C2}" 
        done 
    done
    # iterate over tissue combinations
	#for (( j=$START2; j<$END2; j++ )); do
    #    T2=${TISSUES[$j]}
    #    # do multi    
    #done
done


echo "done"