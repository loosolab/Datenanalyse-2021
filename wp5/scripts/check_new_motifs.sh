#!/bin/bash

# input parameters
DIR="/mnt/workspace_stud/allstud/wp5/runs"
FILE_PATH="/mnt/workspace_stud/allstud/wp5/runs"

# generate and fill check files
for TISSUE in $DIR/*; do
    if [ -d $TISSUE ]; then
        TISSUE=$(echo $TISSUE | rev | cut -d'/' -f-1 | rev)
        # output: tissue_motifs file
        FILE_PATH_TISSUE="${FILE_PATH}/$TISSUE"
        FILE_NAME=check_new_motifs_${TISSUE}.txt
        FILE=${FILE_PATH_TISSUE}/$FILE_NAME
        if ! [ -f $FILE ]; then
            cd $FILE_PATH_TISSUE
            touch $FILE_NAME
            echo "$FILE_NAME was created."
            echo -e "tissue\tcell_type\tsequence\tmotif_name\t" >> $FILE
        fi
    
        TISSUE_AVAILABLE=$(cat $FILE | grep -c "${TISSUE}")
        
        if [ "$TISSUE_AVAILABLE" = 0  ]; then
            for CELL_TYPE in $DIR/$TISSUE/*; do
                # declare the two arrays to store the needed values of the tissue_motifs file
                ID_LIST=()
                CONSENUS_LIST=()
                
                if [ -d $CELL_TYPE/motif_discovery_pipeline ]; then
                    CELL_TYPE=$(echo $CELL_TYPE | rev | cut -d'/' -f-1 | rev)
                    if [ -d $DIR/$TISSUE/$CELL_TYPE/motif_discovery_pipeline/2_discovery ]; then
                        MOTIF_FILE_PATH="$DIR/$TISSUE/$CELL_TYPE/motif_discovery_pipeline/2_discovery"
                        cd $MOTIF_FILE_PATH
                        # variable to ignore the first input line
                        IGNORER=0

                        # read motif_stats_tsv
                        while read line
                        do
                            if [ $IGNORER -eq 1 ]; then
                                # strip the input line and write the values into an array
                                stripped_line=($(echo $line | tr "\t" "\n"))
                                # assign the needed values to variables
                                ID=$(echo ${stripped_line[0]})
                                CONSENUS=$(echo ${stripped_line[11]})
                                # append the variable values to the arrays
                                ID_LIST+=($ID)
                                CONSENUS_LIST+=($CONSENUS)
                            else
                                # start loop after the first input line
                                let IGNORER++
                            fi
                        done < motif_stats.tsv                    
                    fi
                    # counter for indices in arrays
                    COUNTER_ARR=0
                    for elem in "${ID_LIST[@]}"; do
                        echo -e "${TISSUE}\t${CELL_TYPE}\t${ID_LIST[COUNTER_ARR]}\t${CONSENUS_LIST[COUNTER_ARR]}\t" >> $FILE
                        let COUNTER_ARR++
                    done
                fi
            done
        else
            echo "$TISSUE is already analyzed."
        fi
    fi
done


#mnt/workspace_stud/allstud/wp5/runs/liver/cluster20/motif_discovery_pipeline/2_discovery/motif_stats.tsv
#-> TAB seperated
#-> id (erste: 1) = sequence
#-> consensus (letzte: 12) = motif name
#-> ordnen nach Gewebe, Zelltyp, Sequenz, Motivname
#-> Speicherort: für alle Gewebe (runs), für ein Gewebe (tissue)
#-> Barplot (wie oft kommt ein Motiv (Sequenz) vor?)
#-> Dendrogramm (Welche Sequenzen kommen in welchem Gewebe_Zelltyp_Motivname vor)
#-> Analyse: Vergleich ids -> ist etwas gleich? -> wenn ja neue Datei
#-> Aufbau neue Datei:
#>Sequenz1
#    Gewebe1, Zelltyp1, Motivname1
#    Gewebe2, Zelltyp2, Motivname2
#    ...
#>Sequenz2
#    Gewebe1, Zelltyp1, Motivname1
#    Gewebe2, Zelltyp2, Motivname2
#    ...
#-> Barplot (Ranking -> welches Motiv kommt wie oft vor)
#-> Dendrogramm (Welche Sequenzen kommen in welchem Gewebe_Zelltyp_Motivname vor)
