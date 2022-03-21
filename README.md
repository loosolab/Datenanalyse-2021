# Preprocessing of snATAC-seq .FASTQ files
## Introduction
The  provided scripts were written to manage the following steps:
* Create .snap file using [SnapTools](https://github.com/r3fang/SnapTools)
* Using [SnapATAC](https://github.com/r3fang/SnapATAC) for the following steps:
  * Barcode filtering
  * Bin filtering
  * Dimensionality reduction (clustering)
  * Peak calling of each cluster
* Cell type annotation using [python script](https://github.com/loosolab/Datenanalyse-2021/blob/wp2/cell_type_annotation/cell_type_annotation.py), [Uropa](https://github.com/loosolab/UROPA) and [Panglao DB](https://panglaodb.se)

## Setup
Install [SnapTools](https://github.com/r3fang/SnapTools), [SnapATAC](https://github.com/r3fang/SnapATAC) and [Uropa](https://github.com/loosolab/UROPA) as described on the respective homepage.

## IO
**Input:** 
* snATAC-seq .FASTQ file
* .GTF and genome file
* blacklist regions .bed file

**Output:**
* .snap file and .bed file of sample
* plots: barcode filtering, first 50 pairs of eigenvectors, t-sne
* cluster assignment table
* narrow peaks of each cluster
* cell type annotation

## Example
*coming soon*
