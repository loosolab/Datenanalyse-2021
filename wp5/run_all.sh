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
# similarity analysis
${SPATH}/scripts/cluster_all.sh "overall"
echo "Activating plotting environment..."
source /opt/miniconda/bin/activate plotting
${SPATH}/scripts/Eval_Motif_similarity.py --runs-dir $PROJECT_DIR --out "similarity_overall" --motifs "${PROJECT_DIR}/overall_Cluster/motif_comparison_clusters.yml"
echo "Deactivating plotting environment..."
conda deactivate

${SPATH}/scripts/generate_gene_sets.sh
${SPATH}/scripts/generate_gene_sets_TFs.sh
${SPATH}/scripts/compare_gene_sets.py
${SPATH}/scripts/analyse_GO_and_pathway.py

echo "pipeline finished"