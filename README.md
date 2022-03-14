# Datenanalyse-2021

## Introduction
The task of work package 4 was to correlate the chromatin accessibility data from the original ATAC-seq study. The used algorithm was the R package Cicero. 
## Requirements
R Version 4.1.2

## Setup
In order to run Cicero, the packages Monocle3 and Cicero must be installed. 
For installation instructions load the script `Install_MC.R`

## Input Files 

Used input files for this work package are: 

1. `Information about open/closed chromatin`
  expression values in a numeric matrix, rows = peaks, columns = cells
  This file supplies the binary information if a certain peak in a cell is open or closed.
  Provided by WP1, found in directory: `/mnt/workspace_stud/stud2/output/wp4/`
  This file is a sparse matrix with the file extension '_X.mtx'
2. `Cell information`
  data frame, rows = cells (same row names as expression matrix), columns cell attributes 
  This data frame contains the barcodes of the cells. 
  Provided by WP1, found in directory: `/mnt/workspace_stud/stud2/output/wp4/`
  This file is a csv table with the file extension '_obs.csv'
3. `Peak information`
  data frame, rows = peaks, columns = peak attributes
  This data frame supplies information about the open chromatin regions. 
  Provided by WP1, found in directory: `/mnt/workspace_stud/stud2/output/wp4/`
  This file is a csv table with the file extension '_var.csv'
4. `Information about gene attributes`
  The gtf (Gene/General Transfer Format) holds information about gene structure. In this workpackage it is used to connect the peaks to known biological structures: 
  e.g. name of the chromosome, feature type name (e.g. Gene), start postion and end position, strand (+ (forward), - (reverse)), ...  
  direct path to the file: `/mnt/workspace_stud/allstud/homo_sapiens_104.mainChr.gtf`
5. `File with chromosome lengths`
  In order to run `Cicero` the information about the chromosome lengths is required. This information is obtained from the `.fai` file. 
  This data frame holds among others infomration about the name of the sequence (Chromosome) and the total length of if in base pairs, this two specifications are  required.

## Example Case 
