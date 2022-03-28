#!/bin/bash

name_tissue=$1
echo "name_tissue="$name_tissue
#fasta_path="/mnt/data/fastq_data/renlab.sdsc.edu/kai/data/fastq/first/"
# path to files
fasta_path=$4

# fasta files
fasta1="$fasta_path/$name_tissue""_Rep1.demultiplexed"".R1.fastq.gz"
fasta2="$fasta_path/$name_tissue""_Rep1.demultiplexed"".R2.fastq.gz"

echo $fasta1
echo $fasta2
printf "path is :"
echo $3

#exec "bwa"

exec "$3"fasta-to-snap.sh \
--fasta1=$fasta1 \
--fasta2=$fasta2 \
--refgenome=$3genome/hg38.fa \
--refname="hg38" \
--chrom-sizes=$3"genome/hg38.chrom.sizes" \
--short-name=$2 \
--path=$3


