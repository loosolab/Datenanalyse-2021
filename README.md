# Datenanalyse-2021

## Introduction

Workpackage 4 takes a closer look at chromatin peak co-accessibilities. Peaks describe areas of open chromatin. The Cicero algorithm calculates if there are distal peaks to a peak that show the same pattern of open or closed chromatin. Open promoter peaks and corresponding distal peaks are displayed, indicating whether there are open chromatin sites in the distal genome for transcriptional regulation by enhancers or transcription factors.

## Setup
In order to run Cicero, the packages Monocle3 and Cicero must be installed. 
For installation instructions load the script `Install_MC.R`
Once the installation is done, you can run the Script `Run_Cicero.R`, make sure to input the right data if necessary change the input directories. 
__Further information about individual code chunks is provided in the Folder Example Case.__


## Requirements
R Version 4.1.2

## Input Files 

Used input files for this work package are: 

1. #### Information about open/closed chromatin
      expression values in a numeric matrix, rows = peaks, columns = cells.\
      This file supplies the binary information if a certain peak in a cell is open or closed.\
      This file is a sparse matrix with the file extension '_X.mtx' in the R Code the corresponding variable is called 'indata'.\
      Provided by WP1, found in directory: `/mnt/workspace_stud/stud2/output/<tissue_of_interest>/wp4/` 

2. #### Cell information
      data frame, rows = cells (same row names as expression matrix), columns cell attributes.\
      This data frame contains the barcodes of the cells.\
      This file is a csv table with the file extension '_obs.csv' in the R Code the corresponding variable is called 'cellinfo'.\
      Provided by WP1, found in directory: `/mnt/workspace_stud/stud2/output/<tissue_of_interest>/wp4/`

3. #### Peak information
      data frame, rows = peaks, columns = peak attributes.\
      This data frame supplies information about the open chromatin regions.\
      This file is a csv table with the file extension '_var.csv' in the R Code the corresponding variable is called 'peakinfo'.\
      Provided by WP1, found in directory: `/mnt/workspace_stud/stud2/output/<tissue_of_interest>/wp4/`

4. #### Information about gene attributes
      The gtf (Gene/General Transfer Format) holds information about gene structure.\ 
      In this workpackage it is used to connect the peaks to known biological structures:\
      e.g. name of the chromosome, feature type name (e.g. Gene), start postion and end position, strand (+ (forward), - (reverse)), ...\
      direct path to the file: `/mnt/workspace_stud/allstud/homo_sapiens_104.mainChr.gtf`

5. #### File with chromosome lengths
      In order to run `Cicero` the information about the chromosome lengths is required.\
      This information is obtained from the `.fai` file.\
      This data frame holds among others infomration about the name of the sequence (Chromosome) and the total length of if in base pairs, this two specifications are required.

## Output Files

1. In order to connect peaks to a cluster a variable containing a matrix was created called `cl_indata`. This data is the initial input mtx `indata` with rownames of initial peak information file `peakinfo`, but instead of the barcode information of the initial `cellinfo` file this file has the clusters of the cells as column names. 
A function assigns the clusters to the peaks by counting the 1's in the binary `indata.mtx`. This information is stored in the directory `/mnt/workspace_stud/stud10/output/<tissue_of_interest>/cl_indata.mtx`. 

2. `Cicero` provides information about the co-accessibility of chromatin regions. This pairwise comparison is found in the co-accessibility score in the variable  `conns`.  
It's a data frame of pairwise peak comparison with subsequently calculated co-accessibility scores, in a range from -1 to 1. This information is stored in the directory `/mnt/workspace_stud/stud10/output/<tissue_of_interest>/conns.csv`. 

3. The variable `site` stores peaks that are assigned to a promoter. To do so, the initial generated `cds_object` is extended by the information obtained from the .gtf file. The file `cds_sites.csv` is stored in the directory `/mnt/workspace_stud/stud10/output/<tissue_of_interest>/cds_sites.csv`. 

4. `Cicero` also creates gene activity scroes, which combines the score of the regional accessibility with gene expression. This file `cicero_genes_activitiy.mtx` contains genes and the corresponding scores. This information is stored in the directory `/mnt/workspace_stud/stud10/output/<tissue_of_interest>/cicero_gene_activity.mtx`. 


## Example Case 

[Example Run](https://github.com/loosolab/Datenanalyse-2021/blob/wp4/wp4/Example%20Case/Example%20Run.md)
