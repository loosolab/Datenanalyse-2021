#!/bin/sh

# TODO
# check if conda evironment is active
if ["$CONDA_DEFAULT_ENV" != "snakemake"]; then
    echo "Please activate snakemake environment first"
    echo "Abort..."
    exit 1
fi

# remove last slash of directory path if existing
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