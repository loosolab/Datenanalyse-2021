#!/bin/sh

# suffix for found motifs
PATTERN=$1

for FILE in "${@:2}"; do
    sed -E -i "s/motif_([[:digit:]]+)([[:blank:]]+)motif_([[:digit:]]+)/${PATTERN}_\1\2${PATTERN}_\3/" ${FILE}
done