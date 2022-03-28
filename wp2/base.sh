#!/bin/bash

base_path="/mnt/"
path="WP2_workflow_output"
config=""
fastq_path="$base_path""data/fastq_data/renlab.sdsc.edu/kai/data/fastq/first/"
clear=TRUE
tissue=ALL
eigenvector="15"
while [ $# -gt 0 ];do
        case $1 in
              path=*)
                path="${1#*=}"
                ;;
              config=*)
                config="${1#*=}"
                ;;
              clear=*)
                clear="${1#*=}"
                ;;
	      tissue=*)
                tissue="${1#*=}"
                ;;
              eigenvector=*)
                eigenvector="${1#*=}"
                ;;

              *)
                printf "Error: invalid argument: $1"
                exit 1
        esac
        shift
done

help () {

printf "Help TODO \n"


}

#config="$base_path""$path""/wp2_workflow.conf"

create_conf () {
if [ "$config" == "" ];
	then
		printf "config file in $path is created\n"
		config="$base_path""$path""/wp2_workflow.conf"
		touch $config
		echo "#Output directory" > $config
		echo -e "PATH=""$base_path""$path\n" >> $config
		echo "#fastq files path" >> $config
		echo -e "FASTQ_PATH""$fastq_path\n" >> $config
                echo "#number of cores" >> $config
                echo -e "CORES=""5\n" >> $config
                echo "#If TRUE script will process all fastq files" >> $config
                echo -e "ALL=""FALSE\n" >> $config
                echo "#If TRUE script will create .bam and .snap files" >> $config
                echo -e "SNAP=FALSE\n" >> $config
                echo "#If TRUE script will crate Eigenvector Images" >> $config
                echo -e "EIGENVEKTOR=FALSE\n" >> $config
                echo "#IF TRUE script will create Tsne plots" >> $config
                echo -e "TSNE=""FALSE\n" >> $config


	fi

}

print_config () {
echo "#####################"
echo "Create Path: $base_path$path";
echo "config: $config";
echo "fastq files path = $base_path$fastq_path"
echo "#####################"

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

i=0
for f in "${file_arr_unique[@]}"
do

echo "starting for :"$f
IFS="-"
read -a SHORTNAME <<< "$f"
SHORTNAME="${SHORTNAME[-1]}"
SHORTNAME=${SHORTNAME:0:5}${SHORTNAME:12}
echo "$SHORTNAME"
echo "####################################################################################"


printf "testing for file in ""$base_path""WP2_OUTPUT/FINISHED""/""$SHORTNAME""/""$SHORTNAME"".snap"
#snap files
#if [[ ! -f "$base_path""$path""/""$SHORTNAME""/""$SHORTNAME"".snap"  &&  ! -f "$base_path""WP2_OUTPUT/FINISHED""/""$SHORTNAME""/""$SHORTNAME"".snap" ]]; then
if [[ ! -f "$base_path""$path""/""$SHORTNAME""/""$SHORTNAME"".snap"  &&  ! -f "$base_path""WP2_OUTPUT/FINISHED""/""$SHORTNAME""/""$SHORTNAME"".snap" ]]; then
#if [  -f "$base_path""$path""/""$SHORTNAME""/""$SHORTNAME"".snap" ]; then

echo "#####################Start!"$f
#DOSNAP
/mnt/$path/start-fasta-to-snap.sh "$f" "$SHORTNAME" "$base_path$path/"

mkdir "$base_path""$path""/""$SHORTNAME"
mkdir "$base_path""$path""/""$SHORTNAME/plots_and_information"
mkdir "$base_path""$path""/""$SHORTNAME/WP3"
mkdir "$base_path""$path""/""$SHORTNAME/WP6"
mkdir "$base_path""$path""/""$SHORTNAME/annotation"

mv "$base_path""$path""/""$SHORTNAME"".bam" "$base_path""$path""/""$SHORTNAME""/"
#cp "$base_path""$path""/""$SHORTNAME""/""$SHORTNAME"".bam" "$base_path""$path""/""$SHORTNAME""/""$SHORTNAME""_backup"".bam"
mv "$base_path""$path""/""$SHORTNAME"".snap" "$base_path""$path""/""$SHORTNAME""/"
#parsebam
#/mnt/$path/parsebam.py "$base_path""$path""/""$SHORTNAME""/""$SHORTNAME"
mv "$base_path""$path""/""$SHORTNAME""/""$SHORTNAME"".bam" "$base_path""$path""/""$SHORTNAME""/WP3/"
#/mnt/$path/parsebam.py "$base_path""$path""/""$SHORTNAME""/WP3/""$SHORTNAME"
/mnt/$path/parsebam.py "$base_path""$path""/""$SHORTNAME""/WP3/""$SHORTNAME"
mv "$base_path""$path""/""$SHORTNAME""/WP3/""$SHORTNAME""_parsed.bam" "$base_path""$path""/""$SHORTNAME""/WP3/""$SHORTNAME"".bam"


fi
#TODO if condition
#/mnt/$path/parsebam.py "$base_path""$path""/""$SHORTNAME""/WP3/""$SHORTNAME"
#mv "$base_path""$path""/""$SHORTNAME""/WP3/""$SHORTNAME""_parsed.bam" "$base_path""$path""/""$SHORTNAME""/WP3/""$SHORTNAME"".bam"

#eigenvectoren
#if [[ ! -f "$base_path""$path""/""$SHORTNAME""/plots_and_informations/""$SHORTNAME""_Eigenvector_plot.jpeg"  &&  ! -f "$base_path""WP2_OUTPUT/FINISHED""/""$SHORTNAME""/""plots_and_informations/""$SHORTNAME""_Eigenvector_plot.jpeg" ]]; then
#printf "go for eigenvectors""$base_path""$path""/""$SHORTNAME""/plots_and_informations/""$SHORTNAME""_Eigenvector_plot.jpeg"
#Rscript /mnt/$path/preprocessing_WP2.R "$base_path$path/" "$SHORTNAME"
#mv "$base_path""$path""/""$SHORTNAME""_coverage_hist.jpeg" "$base_path""$path""/""$SHORTNAME""/plots_and_informations/""$SHORTNAME""_coverage_hist.jpeg"
#mv "$base_path""$path""/""$SHORTNAME""_Eigenvector_plot.jpeg" "$base_path""$path""/""$SHORTNAME""/plots_and_informations/""$SHORTNAME""_Eigenvector_plot.jpeg"

#fi
#t-sne
if [[ ! -f "$base_path""$path""/""$SHORTNAME""/plots_and_informations/t-sne.png"  &&  ! -f "$base_path""WP2_OUTPUT/FINISHED""/""$SHORTNAME""/plots_and_information/""t-sne.png" ]] ; then
#if [  -f "$base_path""$path""/""$SHORTNAME""/""$SHORTNAME""_t-sne_plot_""$eigencevtor"".jpeg" ]; then
#printf "go for tsne""$base_path""$path""/""$SHORTNAME""/plots_and_informations/""$SHORTNAME""_t-sne_plot_""$eigenvector"".jpeg"
#/mnt/$path/parsebam.py "$base_path""$path""/""$SHORTNAME""/WP3/""$SHORTNAME"
#mv "$base_path""$path""/""$SHORTNAME""/WP3/""$SHORTNAME""_parsed.bam" "$base_path""$path""/""$SHORTNAME""/WP3/""$SHORTNAME"".bam"


Rscript /mnt/$path/preprocessing_WP2_2.R "$base_path$path/" "$SHORTNAME" "$eigenvector"
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



    if [[ "$i" -eq 8 ]]; then
         break
    fi
    ((i++))


done

}



single_workflow () {


echo "starting for :"$tissue

/mnt/testing/start-fasta-to-snap.sh "$tissue"
if [ ! -f "$base_path""$path""/""$tissue""/""$tissue"".snap" ]; then
echo "#####################Start!"$tissue
mkdir "$base_path""$path""/""$tissue"
mkdir "$base_path""$path""/""$tissue/plots_and_informations"
mkdir "$base_path""$path""/""$tissue/WP3"
mkdir "$base_path""$path""/""$tissue/WP6"




mv "$base_path""$path""/""$tissue"".bam" "$base_path""$path""/""$tissue""/"
mv "$base_path""$path""/""$tissue"".snap" "$base_path""$path""/""$tissue""/"
fi
 
#start Rscript for Eigenvectors
#Rscript /mnt/testing/preprocessing_WP2.R "$base_path$path/" "$tissue"



}





#################

#cat /mnt/docker_files/workflow_docker/base.sh
#cat /mnt/docker_files/workflow_docker/base.sh

######
#cat /fasta-to-snap.sh
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


help
create_conf
#echo "$config"
#if  test -f "$config" ; then
#    echo "$config exists."
#    print_config
#    loop_workflow
#else
#   echo "create config"
#   create_conf
#fi

#print_config
loop_workflow

if [ $clear = "TRUE" ];then
	rm "$base_path""$path""/fasta-to-snap.sh"
        rm "$base_path""$path""/start-fasta-to-snap.sh"
        rm "$base_path""$path""/preprocessing_WP2.R"
        rm -r "$base_path""$path""/tmp"

fi

echo "Start..."

#exec "$base_path""$path"/start-fasta-to-snap.sh

#exec which python3
#exec snaptools
#exec /mnt/testing/start-fasta-to-snap.sh
#echo "remove test folder"
#rm -r /mnt/testing
#######################
echo "start r script"
#exec python3 -V
#exec grep
#exec Rscript /mnt/testing/preprocessing_WP2.R "$base_path$path/"

echo "starting...."
#echo "ls"
#find / ENC-1LGRB-194-SM-A9HOW_snATAC_colon_transverse.snap
#echo "ls"
#single_workflow
