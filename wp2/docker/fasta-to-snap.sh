#!/bin/bash

while [ $# -gt 0 ];do
	case $1 in
	      --fasta1=*)
		fasta1="${1#*=}"
		;;
	      --fasta2=*)
		fasta2="${1#*=}"
		;;
	      --refgenome=*)
		refgenome="${1#*=}"
		;;
	      --refname=*)
		refname="${1#*=}"
		;;
	      --chrom-sizes=*)
		chrom_sizes="${1#*=}"
		;;
              --short-name=*)
                short_name="${1#*=}"
                ;;
              --path=*)
                path="${1#*=}"
                ;;


	      --help)
		printf "TODO Usage\n"
		exit 0
		;;
	      *)
		printf "Error: invalid argument: $1"
		exit 1
	esac
	shift
done

### Output definition
out_path="snap_out/"


### Step 1. Index reference genome
bwa=$(which bwa)
bwapath=$(dirname "$bwa")
bwapath=$bwapath"/"
echo $bwapath
#refgenome=refgenome
#refname=refname
name_tmp=$(basename "$fasta1")
name="${name_tmp%.*}"
name="${name%.*}"
name="${name%.*}"
name="${name%.*}"
t_name="${name_tmp%_*}"

printf "bam name is:$t_name"
log="fasta-to-snap.log"

SCRIPTPATH=$(dirname $(realpath -s $0))

printf "LOG FASTA TO SNAP\n"
printf "bwa is = $bwa\n"
printf "refgenome is = $refgenome\n"
printf "refname is = $refname\n\n"
printf "chrom.sizes is $chrom_sizes\n\n"
printf "Fasta1 is $fasta1"
printf "Fasta2 is $fasta2"


index_refgenome () {
printf "starting indexing reference genome...\n"
printf "bwa is = $bwa\n"
printf "refgenome is = $refgenome\n"
printf "refname is = $refname\n\n"
printf "chrom.sizes is $chrom_sizes\n\n"

snaptools index-genome \
	--input-fasta=$refgenome \
	--output-prefix=$refname \
	--aligner=bwa \
	--path-to-aligner=$bwapath \
	--num-threads=5

printf "finished indexing reference genome.\n\n"
}
alignment () {
### Step 2. Alignment

#refgenome=refgenome
#fasta1=fasta1
#fasta2=fasta2

printf "starting alignment...\n"
printf "refgenome is $refgenome\n"
printf "fasta 1 is $fasta1\n"
printf "fasta 2 is $fasta2\n"
printf "output path is : "$path""$short_name
snaptools align-paired-end  \
	--input-reference=$refgenome  \
	--input-fastq1=$fasta1  \
	--input-fastq2=$fasta2  \
	--output-bam=$path$short_name".bam"  \
	--aligner=bwa  \
	--path-to-aligner=$bwapath  \
	--read-fastq-command=zcat  \
	--min-cov=0  \
	--num-threads=5  \
	--if-sort=True  \
	--tmp-folder="$path""tmp"  \
	--overwrite=TRUE

}

pre_processing () {

snaptools snap-pre  \
	--input-file=$path$short_name".bam"  \
	--output-snap=$path$short_name".snap"  \
	--genome-name=hg38  \
	--genome-size=$chrom_sizes  \
	--min-mapq=30  \
	--min-flen=0  \
	--max-flen=1000  \
	--keep-chrm=TRUE  \
	--keep-single=FALSE  \
	--keep-secondary=FALSE  \
	--overwrite=True  \
	--min-cov=100  \
	--verbose=True


}

add_bmat () {

snaptools snap-add-bmat	\
	--snap-file="$path$short_name.snap" 	\
	--bin-size-list 1000 5000 10000	\
	--verbose=True

}

printf "Start indexing ref_genome..." >>$log 
#index_refgenome
printf "finished indexing ref_genome..." >>$log
printf "Start alignment..." >>$log
alignment
echo "finished alignment"
printf "finished alignment..." >>$log
printf "start pre processing..." >>$log
pre_processing
printf "finished pre processing..." >>$log
printf "start add bmat..." >>$log
add_bmat
printf "finished..." >>$log
