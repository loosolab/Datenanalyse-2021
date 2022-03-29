#!/bin/bash

base_path="/mnt/"
path="WP2_workflow_output"
#fastq_path="$base_path""data/fastq_data/renlab.sdsc.edu/kai/data/fastq/first/"
fastq_path="base_path"
clear=TRUE
eigenvector="15"
while [ $# -gt 0 ];do
        case $1 in
              path=*)
                path="${1#*=}"
                ;;
              clear=*)
                clear="${1#*=}"
                ;;
              eigenvector=*)
                eigenvector="${1#*=}"
                ;;
	      fastq_path=*)
                fastq_path="$base_path""${1#*=}"
                ;;
              refgenome=*)
		refgenome="${1#*=}"
		;;
	      chrom_sizes=*)
		chrom_sizes="${1#*=}"
		;;
	      gencode_bed=*)
		gencode_bed="$base_path""${1#*=}"
		;;
	      blacklist=*)
		blacklist="$base_path""${1#*=}"
		;;
	      gtf=*)
		gtf="$base_path""${1#*=}"
		;;
              *)
                printf "Error: invalid argument: $1"
                exit 1
        esac
        shift
done


print_conf () {


		printf "starting with following values \n"
		echo "#Output directory \n"
		echo -e "PATH=""$base_path""$path\n"
		echo "#fastq files path \n"
		echo -e "FASTQ_PATH""$fastq_path \n"
                echo "#number of cores \n"
                echo -e "CORES=""5 \n"
                echo "#Default EIGENVEKTOR value \n"
                echo -e "EIGENVEKTOR=$eigenvector \n"

}


loop_workflow () {

file_arr=()
target="$fastq_path"
echo "$target"
let count=0
for f in "$target"/*
do
    echo $(basename $f)
    name=$(basename $f)
    tissue_name=${name%_Rep*}
    file_arr+=("$tissue_name")
    let count=count+1
done
echo ""
echo "Files in fastq path: $count"

file_arr_unique=( $(printf '%s\n' "${file_arr[@]}" | sort -u) )

for f in "${file_arr_unique[@]}"
do

echo "starting for :"$f
IFS="-"
read -a SHORTNAME <<< "$f"
SHORTNAME="${SHORTNAME[-1]}"
SHORTNAME=${SHORTNAME:0:5}${SHORTNAME:12}
echo "$SHORTNAME"
echo "####################################################################################"


printf "testing for file in ""$base_path""$path""/""$SHORTNAME""/""$SHORTNAME"".snap"
#snap files
if [[ ! -f "$base_path""$path""/""$SHORTNAME""/""$SHORTNAME"".snap"  ]]; then

echo "#####################Start! alignment "$f
#DOSNAP
/mnt/$path/start-fasta-to-snap.sh "$f" "$SHORTNAME" "$base_path$path/" "$fastq_path" "$base_path$refgenome" "$base_path$chrom_sizes"

mkdir "$base_path""$path""/""$SHORTNAME"
mkdir "$base_path""$path""/""$SHORTNAME/plots_and_information"
mkdir "$base_path""$path""/""$SHORTNAME/WP3"
mkdir "$base_path""$path""/""$SHORTNAME/WP6"
mkdir "$base_path""$path""/""$SHORTNAME/annotation"

mv "$base_path""$path""/""$SHORTNAME"".bam" "$base_path""$path""/""$SHORTNAME""/"
#cp "$base_path""$path""/""$SHORTNAME""/""$SHORTNAME"".bam" "$base_path""$path""/""$SHORTNAME""/""$SHORTNAME""_backup"".bam"
mv "$base_path""$path""/""$SHORTNAME"".snap" "$base_path""$path""/""$SHORTNAME""/"
#parsebam
mv "$base_path""$path""/""$SHORTNAME""/""$SHORTNAME"".bam" "$base_path""$path""/""$SHORTNAME""/WP3/"
/mnt/$path/parsebam.py "$base_path""$path""/""$SHORTNAME""/WP3/""$SHORTNAME"
mv "$base_path""$path""/""$SHORTNAME""/WP3/""$SHORTNAME""_parsed.bam" "$base_path""$path""/""$SHORTNAME""/WP3/""$SHORTNAME"".bam"


fi

#t-sne
if [[ ! -f "$base_path""$path""/""$SHORTNAME""/plots_and_informations/t-sne.png" ]] ; then


Rscript /mnt/$path/preprocessing_WP2_2.R "$base_path$path/" "$SHORTNAME" "$eigenvector" "$gencode" "blacklist" "gtf"
mv "$base_path""$path""/t-sne.png" "$base_path""$path""/""$SHORTNAME""/plots_and_information/"
mv "$base_path""$path""/barcodes.png" "$base_path""$path""/""$SHORTNAME""/plots_and_information/"
mv "$base_path""$path""/""$SHORTNAME""_cluster_assignment_table.txt" "$base_path""$path""/""$SHORTNAME""/WP3/"
mv "$base_path""$path""/""$SHORTNAME""peaks.combined.bed" "$base_path""$path""/""$SHORTNAME""/plots_and_information/"

/mnt/$path/narrowpeaks_to_bed.py "$base_path$path/"
find "$base_path""$path""/" -maxdepth 1 -name '*peaks.bed' -exec mv {} "$base_path""$path""/""$SHORTNAME""/WP6/" \;
find "$base_path""$path""/" -maxdepth 1 -name '*.narrowPeak' -exec mv {} "$base_path""$path""/""$SHORTNAME""/plots_and_information/" \;
find "$base_path""$path""/" -maxdepth 1 -name '*.bdg' -exec rm {}  \;

echo "Parse Done"


fi

done

}





echo "create tmp folder..."
mkdir "$base_path"$path
echo  "#logfile" >> "$base_path""$path""/wp2_workflow.log"
mkdir "$base_path"$path"/tmp"
echo "copy scripts..."
cp /fasta-to-snap.sh "$base_path""$path"
cp /start-fasta-to-snap.sh "$base_path""$path"
cp /preprocessing_WP2.R "$base_path""$path"
cp /preprocessing_WP2_2.R "$base_path""$path"
cp /parsebam.py "$base_path""$path"
cp /narrowpeaks_to_bed.py "$base_path""$path"



print_conf
loop_workflow


if [ $clear = "TRUE" ];then
	rm "$base_path""$path""/fasta-to-snap.sh"
        rm "$base_path""$path""/start-fasta-to-snap.sh"
        rm "$base_path""$path""/preprocessing_WP2.R"
        rm -r "$base_path""$path""/tmp"
	rm "$base_path""$path""/parsebam.py"
	rm "$base_path""$path""/narrowpeaks_to_bed.py"
	rm "$base_path""$path""/preprocessing_WP2.R"
	rm "$base_path""$path""/preprocessing_WP2.R"
	rm "$base_path""$path""/wp2_workflow.log"
fi

echo "Finished..."


