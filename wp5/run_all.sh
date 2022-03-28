#!/bin/bash
set -e
## Script to run whole WP5 pipeline

#check for flag:
if ! options=$(getopt -o c -l check-logs -- "$@"  2>/dev/null)
then
    # not correct flag
    echo "Error: script usage: ./$(basename ${0}) [-c | --check-logs]" >&2
    exit 1
fi

set -- $options

while [ $# -gt 0 ]
do
    case $1 in
    -c|--check-logs) check_logs=true ;;
    (--) shift; break;;
    (-?) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    (*) break;;
    esac
    shift
done

# get script path
SPATH=$(dirname $0)
#get config
CONF="${SPATH}/global_vars.cnf"
# clean config 
sed -i 's/\r$//g' $CONF
# read in config
while read LINE; do declare "$LINE"; done < $CONF

# run pipeline steps:
echo "starting pipeline"
${SPATH}/scripts/preparations.sh
${SPATH}/scripts/run_pipeline.sh
if [ "$check_logs" = true ] ; then
    ${SPATH}/scripts/check_logs.sh
fi 
${SPATH}/scripts/renameAllMotifs.sh

# similarity analysis
${SPATH}/scripts/cluster_all.sh "overall"
echo "Activating plotting environment..."
source /opt/miniconda/bin/activate plotting
${SPATH}/scripts/eval_Motif_similarity.py --runs-dir $PROJECT_DIR --out "similarity_overall" --motifs "${PROJECT_DIR}/overall_Cluster/motif_comparison_clusters.yml" --annotation-dir $DATA_PREP_DIR
echo "Deactivating plotting environment..."
conda deactivate

# gene set analysis
if [ $ANN_CHECKER = "yes" ]; then
    ${SPATH}/scripts/generate_gene_sets.sh
    ${SPATH}/scripts/generate_gene_sets_TFs.sh
    ${SPATH}/scripts/compare_gene_sets.py
fi

echo "pipeline finished"