#!/bin/sh

# suffix for found motifs
PATTERN=$1
# motif file
MTF=$2
sed -E -i "s/motif_([[:digit:]]+)       motif_([[:digit:]]+)/motif_\1_${VAR1}    motif_\2_${VAR1}/" ${MTF}
