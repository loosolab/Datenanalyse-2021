# Preprocessing of snATAC-seq FASTQ files
## Introduction
The  provided scripts were written to manage the following steps:
* Create .snap file using [SnapTools](https://github.com/r3fang/SnapTools)
* Using [SnapATAC](https://github.com/r3fang/SnapATAC) for the following steps:
  * Barcode filtering
  * Bin filtering
  * Dimensionality reduction (clustering)
  * Peak calling of each cluster
* Cell type annotation using python script and [Panglao DB](https://panglaodb.se)
