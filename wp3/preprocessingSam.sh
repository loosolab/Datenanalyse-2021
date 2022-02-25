#!/bin/bash
if [ $# = 1 ]
then
    #setting up variables
    inputname=$1
    #replacing "." with " " then creating an array
    inputnamearray=($(echo $inputname | tr "." " "))
    filename=${inputnamearray[0]}
    filenameBam="$filename.bam"
    sortedbam="_sorted.bam"
    filenameSortedBam="$filename$sortedbam"
    samtools view -S -b $inputname > $filenameBam
    samtools sort $filenameBam -o $filenameSortedBam
    samtools index $filenameSortedBam
    rm $filenameBam
    rm $inputname
else
    echo Error: The script expects the filname as argument.
fi
