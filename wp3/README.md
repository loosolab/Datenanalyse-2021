Work Package 3
=======================================

Introduction 
------------

WP3 Contains a script to separate one .bam file of a tissue into seperate files using a cluster assignment table as instructions. It also contains a customised .yaml file for the TOBIAS Snakemake Pipeline, that was used to further process the data from these cluster .bam files.
Afterwards, the footprinting scores are extracted and transcription factor profiles are visualized. 

Installation
------------
Simply clone this repository. Before running the clustering script you will need to also install [pysam](https://pypi.org/project/pysam/) in your python environment. Also, for processing the data, the use of [SAM-Tools](https://www.htslib.org/download/) is recommended.
The post TOBIAS scripts require the python modules argparse, numpy, pandas, scipy, os and matplotlib to function. 

You also need to download the [TOBIAS snakemake repository](https://github.molgen.mpg.de/loosolab/TOBIAS_snakemake) and set up the anaconda environment with the yaml file provided in it like this:

```bash
conda env create -f environments/snakemake.yaml
```
For visualization purposes, you can use [WIlsON](https://academic.oup.com/bioinformatics/article/35/6/1055/5078467). A docker container as well as instructions how to use it are provided [here](https://hub.docker.com/r/loosolab/wilson/).

Workflow description and example
--------------------------------

Converting the .sam file from WP2 to a .bam file and sorting it with samtools. You can use the utility script "preprocessingSam.sh" like this:
```bash
./preprocessingSam.sh sample.sam
```

Alternatively you can do it yourself using these commands:
```bash
samtools view -S -b sample.sam > sample.bam
samtools sort sample.bam -o sample_sorted.bam
samtools index sample_sorted.bam
```
For unsorted .bam files you can alternatively use "preprocessingBam.sh" or leave out the first of the 3 commands above.

Running the clustering script:
```bash
./cluster.py -b input/sample_sorted.bam -t input/cluster_assignment_table.tsv -o clusterBams/
```

The clustering script produces a seperate .bam file for each cluster in the assignment table. It also generates a "snakemakeIn.txt" text file containing the paths and names of the newly created cluster .bam files.

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
After running the snakemake pipeline, further analysis is done with the bindetect_results.txt file for each cluster as well as the bindetect_distances file. These are found in /output/TFBS/. The latter is used to cluster transcription factors into families based on binding motif similarity. This is done via the "TFClustering.py" script. 
```bash
python TFClustering.py bindetect_distances.txt
```
The file "TF_families.tsv" generated by this operation can be used to lower the resolution of the results of the analysis into a more general overview. This is useful because TOBIAS BINDetect is unable to distinguish which transcription factor exactly occupied a position on the DNA based on the provided data, as the scores can only be deduced by the binding motif present on the DNA and thus multiple transcription factors with very similar motifs will also have a similar score. Grouping these transcription factors (which often times have similar functions as well) usually allows for a better overview over the transcription factor spectrum that defines the analysed cluster. 

The binding scores are extracted and brought into a format that allows futher processing via the script "CompareClusterScore.py". 
```bash
python CompareClusterScore.py bindetect_results.txt -f [TF family file]
```
It takes the binding scores for each cluster of the tissue and each transcription factor and writes them down in CLARION file format, which is necessary for heatmap visualization using WIlsON. 
The heatmap gives a visual overview over the data, depicting a scoring profile for each cluster, and shows which transcription factors scored remarkably in some clusters in comparison to the others. After setting up the docker container and opening the localhost in your browser, follow these steps:

Firstly, select "feature" on the top left corner.

![featuure](https://user-images.githubusercontent.com/81377794/160215850-99368edd-89e5-4875-9091-c24488ab6dba.png)

Select or upload the file you want to visualize.

![file](https://user-images.githubusercontent.com/81377794/160215878-8a994827-b225-43e1-88da-95f7e3269b21.png)

Go to "Heatmap > interactive" on the top of the screen.

![heatmap1](https://user-images.githubusercontent.com/81377794/160215931-8d880e31-76f2-4f91-8f88-7af6d66055f7.png)

Select all columns you want to compare on the bottom left, select other desired properties (we reversed the color scheme) nad press "plot" at the bottom.

![columns](https://user-images.githubusercontent.com/81377794/160215976-142c6a55-6211-44e7-9123-6912849effda.png)

You can download the plot at the top right of it.


To gain a better understanding of the transcriptions factors that were remarkably over-represented in certain clusters, the script "DefiningTF.py" extracts the top transcription factors by z score and writes them done in a .tsv file. 
```bash
python DefiningTF.py [CLARION file]
```

This file (as well as those of other tissues at the same time, if that is desired) can be visualized to find out whether transcription factors are defining for multiple cluster or are unique to a certain cell type. This can be done using the script "ProfileSimilarity.py", which generates a bubble plot. The color of the bubble shows in how many clusters this transcription factor made the list of defining transcriptions factors. The bubble size shows how many members are in an examined family (this is done to make it perceptible if a transcription factor family perhaps only entered the defining transcription factors of multiple clusters because it is quite large and different transcription factors in this family might be the reason for this classification). 
```bash
python ProfileSimilarity.py [DefiningTF file] -f [TF family file]
```
The commands shown here are used run the scripts in their default setting for our analysis. They can also all be executed in one command by executing the script "oneTissueAnalysis.sh" while in the diretory of the tissue that contains the snakemake pipeline output (it needs to see the /output/ directory). Giving a name to name the files uniquely according to the tissue that was analyzed is optional, but highly recommended to avoid confusion when working with multiple tissues (escpecially if a comparison between different tissues will be done afterwards).
```bash
bash PATH/TO/FILE/oneTissueAnalysis.sh [TISSUENAME]
```
Different options and flags can be found in the [wiki](https://github.com/loosolab/Datenanalyse-2021/wiki/WP3#wp3-scripts).

The output of the comparison of the clusters in one tissue can also be utilized with these scripts again to compare different tissues to each other.

Known issues and limitations
----------------------------
pysam does only work correctly under Linux operating systems, as it is a requirement for the clustering script it will probably not work on windows.
ProfileSimilarity.py is currently unable to use family clusering when processing multiple input files and the colormap shows decimals when there is little overlap despite it making no biological sense. 

# Example plots

![compareallheatmap](https://user-images.githubusercontent.com/81377794/160215631-d9ecdded-77f4-4875-a37c-547f7958f721.png)

![all](https://user-images.githubusercontent.com/81377794/160215644-570d3df9-996e-4d94-bca0-eb858ecd3274.png)

![liver_heatmap](https://user-images.githubusercontent.com/81377794/160215657-d3d3a5dc-921d-430f-88eb-400076766f46.png)

![liver](https://user-images.githubusercontent.com/81377794/160215663-d39dafaa-e149-4046-b04f-e0f4c13ca299.png)