#!/bin/bash
set -e

## Script to run whole WP5 pipeline
# get script path
SPATH=$(dirname $0)
#get config
CONF="${SPATH}/global_vars.cnf"
# clean config 
sed -i 's/\r$//g' $CONF
# read in config
while read LINE; do declare "$LINE"; done < $CONF

echo "starting pipeline"
# run pipeline steps:
${SPATH}/scripts/generate_configs.sh
${SPATH}/scripts/run_pipeline.sh
${SPATH}/scripts/check_logs.sh
${SPATH}/scripts/renameAllMotifs.sh
${SPATH}/scripts/cluster_all.sh # TODO
${SPATH}/scripts/Eval_Motif_similarity.py # TODO
${SPATH}/scripts/generate_gene_sets.sh
${SPATH}/scripts/generate_gene_sets_TFs.sh
${SPATH}/scripts/compare_gene_sets.py
${SPATH}/scripts/analyse_GO_and_pathway.py

echo "pipeline finished"