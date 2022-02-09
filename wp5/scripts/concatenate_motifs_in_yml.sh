#!/bin/bash

# files
MOTIF_FILE=$1
YML_FILE=$2

# file paths
PATH_MOTIF_FILE="/mnt/workspace_stud/stud12/human_real_data_runs/${MOTIF_FILE}"
PATH_YML_FILE="/mnt/workspace_stud/stud12/human_real_data_runs/${YML_FILE}"

# count the motifs in the motif file
NUMBER_OF_MOTIFS=$(grep -c ^\> $PATH_MOTIF_FILE)

# count the lines in yml file before adding the motifs
LINES_WITHOUT_MOTIFS=$(wc -l < $PATH_YML_FILE)

# read motif file per line, manipulate the motif names and add them at the right position of the yml file
while read LINE; do
    if echo "$LINE" | grep ^\>.*; then
        MODIFIED_LINE=$(echo "$LINE" | sed -e s/\>//)
        sed -i 's,^    motifs:.*$,'"    motifs:\n      - \"$MODIFIED_LINE\""',' $PATH_YML_FILE
    fi
done < $PATH_MOTIF_FILE

# count the lines in yml file after adding the motifs
LINES_WITH_MOTIFS=$(wc -l < $PATH_YML_FILE)
    
# check if the motifs are inserted into yml file
if [ $LINES_WITH_MOTIFS > $LINES_WITHOUT_MOTIFS ] ; then
    if [[ $(( $NUMBER_OF_MOTIFS + $LINES_WITHOUT_MOTIFS )) == $LINES_WITH_MOTIFS ]] ; then
        echo "motifs are added successfully"
    else
        echo "motifs aren't added"
    fi
else
    echo "test failed"
fi
