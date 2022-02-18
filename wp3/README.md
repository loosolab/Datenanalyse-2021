Work Package 3
=======================================

Introduction 
------------

WP3 Contains a script to speperate one .bam file into seperate files using a cluster assignment table as instructions. It also contains a customised .yaml file for the TOBIAS Snakemake Pipeline, that was used to further process the data from these clust .bam files.

Installation
------------
Simply clone this repository. Before running the clustering script you will need to also install [pysam](https://pypi.org/project/pysam/) in your python environment. Aslo for processing the data the use of [SAM-Tools](https://www.htslib.org/download/) is recommended.

You also need to download the [TOBIAS snakemake repository](https://github.molgen.mpg.de/loosolab/TOBIAS_snakemake) and set up the anaconda environment with the yml file provided in it like this:

```bash
conda env create -f environments/snakemake.yaml
```

Workflow description
--------------------

Converting the .sam file from WP2 to a .bam file and sorting it with samtools:
```bash
samtools view -S -b sample.sam > sample.bam
samtools sort sample.bam -o sample_sorted.bam
samtools index sample_sorted.bam
```

Running the clustering script:
```bash
./cluster.py -b ./sample_sorted.bam -t ./cluster_assignment_table.tsv -o ./clusterBams/
```

The clustering script produces a seperate .bam file for each cluster in the assignment table. It also generates a "snakemakeIn.txt" text file containing the paths and names of the cluster .bam files.

You can use the config file provided in this repository (config.yaml).
Paste the text from snakemakeIn.txt into the appropriate section inside the config.yaml and make sure all paths in it lead to the correct files.
Some flags are disabled inside it for this project if you need them you can set them to True instead.

Running the TOBIAS snakemake pipeline:

If it is your first run you might need to update the snakemake version inside your TOBIAS snakemake environment:
```bash
conda activate tobias_snakemake_env
pip install --upgrade snakemake
```
Then you can run the snakemake pipeline:
```bash
conda activate tobias_snakemake_env
snakemake --configfile config.yaml --use-conda --cores [NUMBER_OF_CORES] --conda-prefix /tmp --keep-going
```

Known issues and limitations
----------------------------
pysam does only work correctly under Linux operating systems, as it is a requirement for the clustering script it will probably not work on windows.
