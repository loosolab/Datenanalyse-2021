#!/bin/bash

## Script to start the Motif discovery pipeline for all found config files.
# get script path
SPATH=$(dirname $0)
# read in config
CONF="${SPATH}/../test.conf"
while read LINE; do declare "$LINE"; done < $CONF

# DIR = Directory where the configs are stored
DIR="${PROJECT_DIR}/configs"

# check if conda evironment is active
if  [ ! "$CONDA_DEFAULT_ENV" == "snakemake" ]; then
  echo "Activating snakemake environment..."
  source /opt/miniconda/bin/activate snakemake
fi

# Start pipeline for each config
for CONFIG in $DIR/config_*.yml; do
    dt=$(date '+%d/%m/%Y %H:%M:%S');
        echo "starting ${CONFIG} at ${dt}"
    {   
       snakemake --configfile $CONFIG --cores 5 --use-conda --conda-frontend mamba 

    } || {
        # catch error
        echo "Error with ${CONFIG}"
    }
done

echo "Done"
