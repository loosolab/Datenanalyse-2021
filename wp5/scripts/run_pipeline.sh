#!/bin/bash

## Script to start the Motif discovery pipeline for all found config files.

# check if conda evironment is active
if  [ ! "$CONDA_DEFAULT_ENV" == "snakemake" ]; then
  echo "Activating snakemake environment..."
  source /opt/miniconda/bin/activate snakemake
fi

# remove last slash of directory path if existing
# DIR = Directory where the configs are stored
if [[ $1 == */ ]]; then
  DIR=${1::-1};
else
  DIR=$1; 
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
