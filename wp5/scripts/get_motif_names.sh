#!/bin/bash

FILE=$1
OUT=$2
FORMAT=$3

if [ $FORMAT == "MEME" ] ; then
    grep "MOTIF" $FILE | cut -f2 | awk '{print ">"$0}'> $OUT
else
    grep ">" $FILE | cut -f1 > $OUT
fi