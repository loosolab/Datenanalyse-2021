#!/bin/bash

#script needs to be executed in the tissue directory, so that you can see the output/ directory of the snakemake pipeline
#for our analysis, this is /mnt/workplace_stud/allstud/wp3/new_tissiues/[TISSUENAME]/

#input variable for naming the output, giving it a name to identify the tissue that was analyzed
NAME=$1

#names of the output files
CCNAME="_ClusterComparison.clarion"
DTFNAME="_DefiningTFsPerCluster.tsv"
PSNAME="_Similarity.tsv"

#creating a directory for the analysis files
mkdir Analysis

#executing the analysis scripts
python /mnt/workspace_stud/stud8/Datenanalyse-2021/wp3/TFClustering.py output/TFBS/bindetect_distances.txt

python /mnt/workspace_stud/stud8/Datenanalyse-2021/wp3/CompareClusterScore.py output/TFBS/bindetect_results.txt -f TF_families.tsv

python /mnt/workspace_stud/stud8/Datenanalyse-2021/wp3/DefiningTF.py ClusterComparison.clarion

python //mnt/workspace_stud/stud8/Datenanalyse-2021/wp3/ProfileSimilarity.py importantTFsPerCluster.tsv -f TF_families.tsv

#renamning the output files so they are unique to the tissue
mv ClusterComparison.clarion $NAME$CCNAME
mv importantTFsPerCluster.tsv $NAME$DTFNAME
mv similarity.tsv $NAME$PSNAME

#moving the output files to the analysis diretory
mv $NAME$CCNAME Analysis/
mv $NAME$DTFNAME Analysis/
mv $NAME$PSNAME Analysis/
mv TF_families.tsv Analysis/

