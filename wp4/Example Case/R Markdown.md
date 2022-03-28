# Introduction 

Chromatin, which consists of DNA, is a dynamic structure. Regulation of transcription is based on the interaction between the structure of chromatin and the recruitment of numberous transcription factors, proximal promoter elements and upstream activator sequences. 
Accessbile chromatin is crutial for transcriptional regulation. This accessibility is marked by DNA methylation and histone modification. 
Accessbility can be modified by both pathogenic and environmental factors, indicating position of regulatory regions and by that reflecting regulation of cell behavior.  

To visualize regions of accessible chromatin in real time (genome-wide) the method Assay for Transpoase Accessible Chromatin with high-throughput sequencing (ATAC-seq) was generated[<sup>1</sup>](#fn1). 

`Cicero` was developed as algorithm to create a link between this chromatin accessibility of regulatory elements to their target genes and is by that creating a prediction of gene expression. 

# Looso-Lab Project 

This student project was created by the team around Mario Looso from the MPI in Bad-Nauheim. In the course Biodatenanalyse, data of an already conducted study[<sup>2</sup>](#fn2) should be bioinformatically processed. The team chose the study "A cell atlas of chromatin accessibility across 25 adult human tissues" from Zhang et. al. 
This study was conducted by the usage of 70 bio-samples, with 25 different tissue types, obtained from four donors. The clustering of 472.373 nuclei resulted in 54 obvious cell types. 
The global goal of this student project was to work with this data gained from sci-ATAC-seq assays. 
By splitting the workload into individual work packages, the raw data was further processed, meeting a specific biological question per work package. 

This work package looks at the co-accessible chromatin regions, creating a score to link chromatin regions with same opening scheme to distal regions in the upstream genome and annotating the open regions to known promoters, transcription start sites or transcription factors. The software `Cicero` was used to get to the bottom of this question. 

# Cicero Introduction

Cicero[<sup>3</sup>](#fn3) is an algorithm identifying co-accessibility pairs of DNA elements whilst connecting regulatory elements to their putative target genes based on dynamics of accessibility of linked distal elements. 

Other approaches differ from Cicero which isn't using bulk chromatin accessibility generated over many cell lines and tissues but working with single cell chromatin accessibility data from a single experiment. `Cicero` has to be robust to the sparsity of the data. 

#### Algorithm[<sup>4</sup>](#fn4) 
By sampling and aggregation of similar cells in groups `Cicero` quantifies correlations between putative regulatory elements. Based on this quantification `Cicero` links regulatory elements to their target genes using machine learning. The algorithm is applicable to every organism and any cell types. The challenges of building a genome-wide cis-regulatory map `Cicero` is facing are the following: 

1. raw correlations are driven by technical features (e.g. read depth/cell) 

2. insufficient observations to estimate correlations between billions of pairs of sites. 

3. singe cell ATAC-seq data is very sparse 

4. while the accessibility of distal elements might be correlated with their target promoters, very distant or interchromosal pair of sites are also correlated because they are part of the same regulatory program 

The user provides `Cicero` with clustered cells as input. The algorithm then creates a great amount of cell groups, each group containing 50 cells in similar position in clustering. This ensures the overcome of the sparsity of data, it furthermore aggregates accessible profiles for cells in groups to produce counts that subtract effects of technical variables and it measures correlations in accessibility between all pairs of sites inside a 500 kb frame. The output of `Cicero` consists of this correlations, the co-accessibility-score. 

![Figure3.jpg](attachment:Figure3.jpg)
***
<center>
An overview of the Cicero algorithm
</center>

`Cicero` can also identify CCANs: cis-co-accessibility networks which are modules of sites highly co-accessible with one another with help of a community detection algorithm. Adapted from the definition of a chromatin hub CCANs should meet the following criteria. 

1. They should be located in close physical proximity, closer than expected based on their proximity in the linear genome (further analysed and verified with data from ChIA-PET analysis  (https://www.sciencedirect.com/science/article/abs/pii/S1046202312002204)). 

2. They should interact with common groups of protein complexes. 

3. The epigenetical modifications should occur at similar times. 

4. They should regulate genes with promoters in the hub substantially. 

Links created by `Cicero`are also mediated by interacting transcription factors (TFs). 

#### Conclusion  
In contrast to other approaches `Cicero`operates with single-cell data, therefore avoiding bulk average effects while analyzing. Analyses with `Cicero` can accelerate the quantitative understanding of eukaryotic gene regulation. Additionally it may ease identification of target genes of non-coding variants of genome wide association signals. The algorithm provides an effective resource to generate links between regulatory elements and target genes in tissues or cell type by using data from a single cell experiment. The defined chromatin hubs help the construction of gene expression dynamic models, furthermore the identification  of genes in which dysregulation is subject to genome-wide association.
As this epigenetic field evolves cell atlases defining each cell type and its molecular profile regulatory maps will be essential for understanding the gene expression program in both illness and health. 

#### Limitations 
The limitation of `Cicero`is mainly based on the putative type of it's generated connections. Further experiments are necessary to determine whether a linked distal DNA element is essential for regulatory influence or sufficient. 

***
# Example Run 

## Section 1: Read in data and create CDS object

`Cicero` is a R package for single cell-ATAC-Seq analysis and an extension of the `Monocle` R package. Both created by the Cole Trapnell Lab. For further documentation on `Cicero` visit https://cole-trapnell-lab.github.io/cicero-release/. For installation instruction visit the R Script `Install_MC.R`. 
The code for this section is found in `Create_CDS.R`. 
After installation of both packages, the librarys has to be loaded to use associated functions.


```R
library(Matrix)
library(stringr)
library(monocle3)
library(cicero)
```

    Loading required package: Biobase
    
    Loading required package: BiocGenerics
    
    
    Attaching package: ‘BiocGenerics’
    
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, sd, var, xtabs
    
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, append, as.data.frame, basename, cbind, colnames,
        dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
        grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
        union, unique, unsplit, which.max, which.min
    
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    Loading required package: SingleCellExperiment
    
    Loading required package: SummarizedExperiment
    
    Loading required package: MatrixGenerics
    
    Loading required package: matrixStats
    
    
    Attaching package: ‘matrixStats’
    
    
    The following objects are masked from ‘package:Biobase’:
    
        anyMissing, rowMedians
    
    
    
    Attaching package: ‘MatrixGenerics’
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars
    
    
    The following object is masked from ‘package:Biobase’:
    
        rowMedians
    
    
    Loading required package: GenomicRanges
    
    Loading required package: stats4
    
    Loading required package: S4Vectors
    
    
    Attaching package: ‘S4Vectors’
    
    
    The following objects are masked from ‘package:Matrix’:
    
        expand, unname
    
    
    The following objects are masked from ‘package:base’:
    
        expand.grid, I, unname
    
    
    Loading required package: IRanges
    
    Loading required package: GenomeInfoDb
    
    
    Attaching package: ‘monocle3’
    
    
    The following objects are masked from ‘package:Biobase’:
    
        exprs, fData, fData<-, pData, pData<-
    
    
    Loading required package: Gviz
    
    Loading required package: grid
    


***
The first step that must be taken is to create a object of the ***cell_data_set (CDS)*** class. This format provides an interface defined from the Bioconductor ***SingleCellExperiment*** class. Three input files are required. 

1. ***expression_matrix***:  expression values in a numeric matrix, rows = genes, columns = cells. 

2. ***cell_metadata***: data frame, rows = cells (same row names as expression matrix), columns cell attributes 

3. ***gene_metadata***: data frame, rows = features (e.g. genes), columns = gene attributes

Since this project uses 10X scATAC-Seq data, a few more steps are required. 
These files are extracted from a .h5ad file using `Scanpy`. Extraction instructions from Anndata can be found in the script `Extract files.jyptnm`. However this data was extracted and provided from Workpackage 1. 

This ***expression_matrix*** contains binary information: is this site open in this cell. In other words: have we found an open chromatin site or have we found a fragment in the sc-ATAC-Seq analysis at this site in this cell. 
We read in the matrix data using the Matrix package.
Binarizing the matrix prevents other entries than 0 and 1. There might be situations in which the transpoase cuts 2 times within one big peak. Cicero can just read matrices who are binarized, therefore matrices only with entries 0 and 1. 


```R
indata <- Matrix::readMM("/mnt/workspace_stud/stud2/output/esophagus_mucosa/wp4/esophagus_mucosa_X.mtx")
indata@x[indata@x > 0] <- 1
indata <- t(indata)
# filter out peaks with less than 30 cells 
indata <- indata[rowSums(indata) > 30,]
head(indata)
```


    6 x 28070 sparse Matrix of class "dgTMatrix"
                                                                                   
    [1,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ......
    [2,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ......
    [3,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ......
    [4,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ......
    [5,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ......
    [6,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ......
    
     .....suppressing 28036 columns in show(); maybe adjust 'options(max.print= *, width = *)'
     ..............................


The information the data frame ***cellinfo*** provides is the information of the assignment of a cell into a cluster. Read in the file with the cell-information by using the function `read.csv()`, specifying the correct separator. 


```R
cellinfo <- read.csv("/mnt/workspace_stud/stud2/output/esophagus_mucosa/wp4/esophagus_mucosa_obs.csv", 
                      sep = '+', row.names = NULL)
```

Change the row names of this data frame in to the first column of the data frame with the function `row.names()` and call this row ***cells*** then by using `names() <- "cells`. Now the first row is called ***cells*** and contains the specific assignment (barcode) of the cell to a cluster. This row will be the columns of the ***expression_matrix***. Further information of the ***cellinfo*** data frame are the allocations to a ***cluster*** and the ***tissue***. 


```R
# Changing the format into data.frame to perform strsplit().
cellinfo <- data.frame(cellinfo)

# Split column barcode by character "\t".
cellinfo[c('code','cluster')] <- str_split_fixed(cellinfo$barcode, '\t', 3)

# Split column code by character "-".
cellinfo[c('barcode', 'cl')] <- str_split_fixed(cellinfo$code, '-', 2)

# Change columns and their column names.
cellinfo <- cellinfo[,  c("barcode", "cluster", "row.names")]
names(cellinfo) <- c("barcode", "cluster", "tissue")

# Set rownames of cellinfo data frame.
rownames(cellinfo) <- NULL
suppressWarnings(rownames(cellinfo) <- cellinfo$barcode)
head(cellinfo)
```


<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>barcode</th><th scope=col>cluster</th><th scope=col>tissue</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AACGACGTGTAATGGTTCCCTT</th><td>AACGACGTGTAATGGTTCCCTT</td><td>0 </td><td>esophagus_mucosa_SM-A9HOR_1</td></tr>
	<tr><th scope=row>AAGCAAAGTCTTCGGTTCCCAC</th><td>AAGCAAAGTCTTCGGTTCCCAC</td><td>5 </td><td>esophagus_mucosa_SM-A9HOR_1</td></tr>
	<tr><th scope=row>AAGTCCTTAGACTCTTCCTCAT</th><td>AAGTCCTTAGACTCTTCCTCAT</td><td>0 </td><td>esophagus_mucosa_SM-A9HOR_1</td></tr>
	<tr><th scope=row>AAGTCCTTAGCAGTGATTCTGT</th><td>AAGTCCTTAGCAGTGATTCTGT</td><td>0 </td><td>esophagus_mucosa_SM-A9HOR_1</td></tr>
	<tr><th scope=row>AAGTCCTTAGGTACAACTCTAG</th><td>AAGTCCTTAGGTACAACTCTAG</td><td>10</td><td>esophagus_mucosa_SM-A9HOR_1</td></tr>
	<tr><th scope=row>ACAATGCTCCTTCCGCTGATAT</th><td>ACAATGCTCCTTCCGCTGATAT</td><td>0 </td><td>esophagus_mucosa_SM-A9HOR_1</td></tr>
</tbody>
</table>



In the Section Biological Approaches the information of the barcode combined with the cluster is required. The variable ***bar_cluster*** stores this information as a vector and is stored in the directory `/output/<tissue>/bar.cluster.tsv`.


```R
bar_cluster <- paste(cellinfo$barcode, cellinfo$cluster, sep="_")
write.csv(bar_cluster, '/mnt/workspace_stud/stud10/output/esophagus_mucosa/bar_cluster.tsv')
```

The data frame ***peakinfo*** contains the information about the peaks.
For the purpose of pre-filtering the column ***shared_cells*** is extracted and only cells with a score higher than 30 are taken into account. 


```R
# format peak info
peakinfo <- read.csv("/mnt/workspace_stud/stud2/output/esophagus_mucosa/wp4/esophagus_mucosa_var.csv",
                     header = TRUE, sep = '_', row.names = NULL, col.names= c("chr", "bp1", "shared", "cells.variability", "score"))
```


```R
# Changing the format into data.frame to perform strsplit().
peakinfo <- data.frame(peakinfo)

# Split column shared by character "\t".
peakinfo[c('bp2', 'number')] <- str_split_fixed(peakinfo$shared, '\t', 2)

# Split column number by character "\t".
peakinfo[c('shared_cells', 'number')] <- str_split_fixed(peakinfo$number, '\t', 2)

# In order to filter the peaks, shared_cells-col must be a numeric information and not stored as character. 
peakinfo$shared_cells <- as.numeric(as.character(peakinfo$shared_cells))

# filter out counts (stored in shared_cells) with less than 30 cells.
peakinfo <- peakinfo[peakinfo$shared_cells > 30,]
```

Create a new column with the information which chromosome, where the peak starts and where the peak ends, delimited through the `_` character. This row will then be the rows of the ***expression_matrix***. 


```R
# Creating new column and paste this concatenated string as new row names. 
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")

# Set rownames of the dataframe
row.names(peakinfo) <- peakinfo$site_name

# Discard unused columns. 
peakinfo <- peakinfo[, c("site_name", "chr", "bp1", "bp2")]
head(peakinfo)
```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>site_name</th><th scope=col>chr</th><th scope=col>bp1</th><th scope=col>bp2</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>chr1_181273_181673</th><td>chr1_181273_181673</td><td>chr1</td><td>181273</td><td>181673</td></tr>
	<tr><th scope=row>chr1_184281_184681</th><td>chr1_184281_184681</td><td>chr1</td><td>184281</td><td>184681</td></tr>
	<tr><th scope=row>chr1_267803_268203</th><td>chr1_267803_268203</td><td>chr1</td><td>267803</td><td>268203</td></tr>
	<tr><th scope=row>chr1_629065_629465</th><td>chr1_629065_629465</td><td>chr1</td><td>629065</td><td>629465</td></tr>
	<tr><th scope=row>chr1_632120_632520</th><td>chr1_632120_632520</td><td>chr1</td><td>632120</td><td>632520</td></tr>
	<tr><th scope=row>chr1_779618_780018</th><td>chr1_779618_780018</td><td>chr1</td><td>779618</td><td>780018</td></tr>
</tbody>
</table>



In order to connect peaks with associated clusters a variable ***cl_indata*** is initialized. ***cl_indata*** is ***indata***, with differing colnames: the cluster information of ***cellinfo***. Store this matrix in `/mnt/workspace_stud/stud10` to work with that on another point of time. 

Prepare the ***expression_matrix***: Set the column names by using `colnames()`and the row names using `row.names()` to assign new values created above. 


```R
row.names(indata) <- rownames(peakinfo)
colnames(indata) <- rownames(cellinfo)

cl_indata <- indata
row.names(cl_indata) <- rownames(peakinfo)
colnames(cl_indata) <- bar_cluster

write(x = rownames(cl_indata), file = "/mnt/workspace_stud/stud10/output/esophagus_mucosa/peaks.tsv")
write(x = colnames(cl_indata), file = "/mnt/workspace_stud/stud10/output/esophagus_mucosa/cluster.tsv")

cl_indata <- Matrix(cl_indata , sparse = T )
head(cl_indata)
writeMM(obj = cl_indata, file="/mnt/workspace_stud/stud10/output/esophagus_mucosa/cl_indata_esophagus_mucosa.mtx")
```

       [[ suppressing 34 column names 'AACGACGTGTAATGGTTCCCTT_0', 'AAGCAAAGTCTTCGGTTCCCAC_5', 'AAGTCCTTAGACTCTTCCTCAT_0' ... ]]
    



    6 x 28070 sparse Matrix of class "dgTMatrix"
                                                                                  
    chr1_181273_181673 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    chr1_184281_184681 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    chr1_267803_268203 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    chr1_629065_629465 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    chr1_632120_632520 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    chr1_779618_780018 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                     
    chr1_181273_181673 . . . . ......
    chr1_184281_184681 . . . . ......
    chr1_267803_268203 . . . . ......
    chr1_629065_629465 . . . . ......
    chr1_632120_632520 . . . . ......
    chr1_779618_780018 . . . . ......
    
     .....suppressing 28036 columns in show(); maybe adjust 'options(max.print= *, width = *)'
     ..............................



    NULL


The function `new_cell_data_set()` creates a CDS object with the arguments: expression_data (***indata***), cell_metadata (***cellinfo***) and gene_metadata (***peakinfo***).


```R
cds_object <- suppressWarnings(new_cell_data_set(indata,
                                cell_metadata = cellinfo, 
                                gene_metadata = peakinfo))
```

`detect_genes()` is a Monocle3 function which counts how many cells are expressed and detectable above a minimum threshold. The results of this functions are added in the columns ***num_cell_expressed*** and ***num_genes_expressed*** in the newly created cds_object.


```R
cds_object <- monocle3::detect_genes(cds_object)
```

This chunk of code ensures that there are no peaks included with zero reads. In our project, this is intercepted by the work of the upstream work packages.


```R
cds_object <- cds_object[Matrix::rowSums(exprs(cds_object)) != 0,] 
```

To clear out memory several variables can be deleted. `save.image()` is a R function to save all variables created (and not deleted) so far. `ls()` displays which variables are stored. Then the kernel can be restarted to start off with newly acquired memory.  


```R
# one can discard variables indata, cellinfo, peakinfo  & cl_indata  to free memory 
rm("indata")
rm("cellinfo")
rm("peakinfo")
rm("cl_indata")
rm("bar_cluster")

save.image()
ls()
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'cds_object'</li><li>'cicero_cds'</li></ol>



***
## Section 2: Preprocess CDS and make Cicero CDS

The code to this section is found in `Create_Cicero_CDS.R`. Firstly the librarys have to be loaded again to make use of associated functions. 


```R
library(ggplot2)
library(dplyr)
library(monocle3)
library(cicero)
```

    
    Attaching package: ‘dplyr’
    
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    
    Loading required package: Biobase
    
    Loading required package: BiocGenerics
    
    
    Attaching package: ‘BiocGenerics’
    
    
    The following objects are masked from ‘package:dplyr’:
    
        combine, intersect, setdiff, union
    
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, sd, var, xtabs
    
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, append, as.data.frame, basename, cbind, colnames,
        dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
        grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
        union, unique, unsplit, which.max, which.min
    
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    Loading required package: SingleCellExperiment
    
    Loading required package: SummarizedExperiment
    
    Loading required package: MatrixGenerics
    
    Loading required package: matrixStats
    
    
    Attaching package: ‘matrixStats’
    
    
    The following objects are masked from ‘package:Biobase’:
    
        anyMissing, rowMedians
    
    
    The following object is masked from ‘package:dplyr’:
    
        count
    
    
    
    Attaching package: ‘MatrixGenerics’
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars
    
    
    The following object is masked from ‘package:Biobase’:
    
        rowMedians
    
    
    Loading required package: GenomicRanges
    
    Loading required package: stats4
    
    Loading required package: S4Vectors
    
    
    Attaching package: ‘S4Vectors’
    
    
    The following objects are masked from ‘package:dplyr’:
    
        first, rename
    
    
    The following objects are masked from ‘package:base’:
    
        expand.grid, I, unname
    
    
    Loading required package: IRanges
    
    
    Attaching package: ‘IRanges’
    
    
    The following objects are masked from ‘package:dplyr’:
    
        collapse, desc, slice
    
    
    Loading required package: GenomeInfoDb
    
    
    Attaching package: ‘monocle3’
    
    
    The following objects are masked from ‘package:Biobase’:
    
        exprs, fData, fData<-, pData, pData<-
    
    
    Loading required package: Gviz
    
    Loading required package: grid
    


The following code chunk is the process of trajectory building of Monocle3. 
In response to stimuli, during development and within their life span cells transition from one functional stat to another. In different states they will express different genes, displaying a dynamic repertoire of proteins and metabolites. 
Moving between these states, they go trough a process of transcriptional reconfiguration. The discretion of these different states are modeled by the Monocle algorithm, learning the sequence of gene expression changes each cell must go through as part of this dynamic biological process. Moncole places the cells in the proper position in the trajectory. The more dimensions there are, the harder is the learning of the trajectory. Monocle provides different algorithms to reduce dimensionality. 

The function `detect_genes()` counts for each gene in a CDS object, how many cells are expressed (above a threshold). Furthermore it counts the number of genes for each cell that are detectable. This information is found in the added columns num_genes_expressed in the rowData and colData tables respectively.
The statistical function `estimate_size_factors()` displays a set of statistical data of the CDS object, like class, rownames, colnames, reducedDimNames and many more. 
Normalization and pre-processing steps are important to normalize data by log and size factor in order to reduce depth differences. The arguments of the function performing this `preprocess_cds()` are the CDS upon which this operation is performed, as well as the method. In the case of ATAC-Seq, LSI (latent semantic indexing) is applied and will set the parameter ***reducedDimNames*** of the ***cds_object*** to 1, LSI. 


```R
cds_object <- detect_genes(cds_object)
cds_object <- estimate_size_factors(cds_object)
cds_object <- preprocess_cds(cds_object, 
                             method = "LSI")
```

The function `reduce_dimension()` holds the following arguments: the CDS upon which the operation is performed, the reduction method, specifying the algorithm used for dimensional reduction, the preprocessing method, indicating which method is used on the data, the ***umap.min_dist***, a numeric indication for the minimum distance to be passed to the UMAP function, as well as the ***umap.n_neighbors***, an integer indicating the number of neighbors during the graph construction. Raising the parameter reducedDimNames of the cds_object to 2, LSI and UMAP. 


```R
cds_object <- reduce_dimension(cds_object, 
                               reduction_method = 'UMAP', 
                               preprocess_method = "LSI", 
                               umap.min_dist = 0.4,
                               umap.n_neighbors=15)
```

A UMAP plot arranges data in low-dimensional space by constructing a high dimensional graph representation of data with subsequent optimization of a low-dimensional graph as structurally similar as possible. UMAP builds a weighted graph, with edge weights showing the likelihood that two points are connected. A radius is spreads outward from each point, connecting to points where radii overlap. This radius is chosen based on the distance to the nth nearest neighbor. At the end each point must be connected at least to its closest neighbor, preserving the balance of the global structure[<sup>5</sup>](#fn5). 

The function `cluster_cells()` returns a CDS object with internally stored cluster assignments using Louvain/Leiden community detection. 


```R
cds_object <- cluster_cells(cds_object)
plot_cells(cds_object, 
           color_cells_by = 'cluster',
           group_label_size = 4,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_groups_by_cluster = TRUE)
```


    
![png](output_38_0.png)
    


Because of the sparsity of single-cell chromatin accessibility data, an estimation of the co-accessibility scores requires an aggregation of similar cells in order to create more dense count data. The approach of k-nearest neighbor, creation of overlapping sets of cells, is used by Cicero to do so. Cicero then constructs these sets while using a reduced dimension coordinate map, representing cell similarity. 

The method `reducedDims()` is to get or set reduction results in a SingleCellExperiment object. Each row of this reduced dimension results are expected to correspond to a column of the SingleCellExperiment object. 
The dimension map ***umap_cords*** is a data frame, in which the row names match the cell IDs in the CDS object and the columns are the coordinates of the reduced dimension object. This data frame represents the coordinates of each cell in reduced dimension space (2-3 dimensions).
The input of the function `make_cicero_cds()` is then the input CDS and the previously created map of reduced dimensions (***umap_coords***), the output is an aggregated input CDS for Cicero. The aggregation is done by the usage of the k-nearest-neighbors graph. 


```R
umap_coords <- reducedDims(cds_object)$UMAP
head(umap_coords)
cicero_cds <- make_cicero_cds(cds_object, reduced_coordinates = umap_coords)
```


<table class="dataframe">
<caption>A matrix: 6 × 2 of type dbl</caption>
<tbody>
	<tr><th scope=row>AACGACGTGTAATGGTTCCCTT</th><td> 0.21205975</td><td>10.769368</td></tr>
	<tr><th scope=row>AAGCAAAGTCTTCGGTTCCCAC</th><td> 1.80137184</td><td> 4.796061</td></tr>
	<tr><th scope=row>AAGTCCTTAGACTCTTCCTCAT</th><td> 0.04174618</td><td>10.577109</td></tr>
	<tr><th scope=row>AAGTCCTTAGCAGTGATTCTGT</th><td>-0.15902084</td><td>11.313312</td></tr>
	<tr><th scope=row>AAGTCCTTAGGTACAACTCTAG</th><td>-0.95649895</td><td> 9.369722</td></tr>
	<tr><th scope=row>ACAATGCTCCTTCCGCTGATAT</th><td>-0.75647167</td><td>11.586324</td></tr>
</tbody>
</table>



    Overlap QC metrics:
    Cells per bin: 50
    Maximum shared cells bin-bin: 44
    Mean shared cells bin-bin: 0.0882025563375618
    Median shared cells bin-bin: 0
    


Display the attributes of ***cicero_cds***. 


```R
suppressWarnings(cicero_cds)
```


    class: cell_data_set 
    dim: 261604 4076 
    metadata(1): cds_version
    assays(1): counts
    rownames(261604): chr1_181273_181673 chr1_184281_184681 ...
      chrY_19737748_19738148 chrY_20575864_20576264
    rowData names(5): site_name chr bp1 bp2 num_cells_expressed
    colnames(4076): agg2986 agg14134 ... agg17220 agg616
    colData names(3): agg_cell Size_Factor num_genes_expressed
    reducedDimNames(0):
    mainExpName: NULL
    altExpNames(0):



```R
# one can discard the variable umap_coords
rm("umap_coords")

save.image()
ls()
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'cds_object'</li><li>'cicero_cds'</li></ol>



***
## Section 3: Create outputfiles and run Cicero

After restarting the kernel the librarys have to be loaded in order to use the associated functions. 
This code is also stored in the file `Run_Cicero.R`. 


```R
library(ggplot2)
library(dplyr)
library(monocle3)
library(cicero)
library(seqTools)
library(stringr)
library(Matrix)
```

    
    Attaching package: ‘dplyr’
    
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    
    Loading required package: Biobase
    
    Loading required package: BiocGenerics
    
    
    Attaching package: ‘BiocGenerics’
    
    
    The following objects are masked from ‘package:dplyr’:
    
        combine, intersect, setdiff, union
    
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, sd, var, xtabs
    
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, append, as.data.frame, basename, cbind, colnames,
        dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
        grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
        union, unique, unsplit, which.max, which.min
    
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    Loading required package: SingleCellExperiment
    
    Loading required package: SummarizedExperiment
    
    Loading required package: MatrixGenerics
    
    Loading required package: matrixStats
    
    
    Attaching package: ‘matrixStats’
    
    
    The following objects are masked from ‘package:Biobase’:
    
        anyMissing, rowMedians
    
    
    The following object is masked from ‘package:dplyr’:
    
        count
    
    
    
    Attaching package: ‘MatrixGenerics’
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars
    
    
    The following object is masked from ‘package:Biobase’:
    
        rowMedians
    
    
    Loading required package: GenomicRanges
    
    Loading required package: stats4
    
    Loading required package: S4Vectors
    
    
    Attaching package: ‘S4Vectors’
    
    
    The following objects are masked from ‘package:dplyr’:
    
        first, rename
    
    
    The following objects are masked from ‘package:base’:
    
        expand.grid, I, unname
    
    
    Loading required package: IRanges
    
    
    Attaching package: ‘IRanges’
    
    
    The following objects are masked from ‘package:dplyr’:
    
        collapse, desc, slice
    
    
    Loading required package: GenomeInfoDb
    
    
    Attaching package: ‘monocle3’
    
    
    The following objects are masked from ‘package:Biobase’:
    
        exprs, fData, fData<-, pData, pData<-
    
    
    Loading required package: Gviz
    
    Loading required package: grid
    
    Loading required package: zlibbioc
    
    
    Attaching package: 'Matrix'
    
    
    The following object is masked from 'package:S4Vectors':
    
        expand
    
    


The GTF (General Transfer Format) file is a variation of the GFF (General Feature Format) File Format. It consists of one line per feature, displaying 9 columns of data, such as ***seqname*** = name of chromosome/scaffold, ***feature*** = feature type name (e.g. Gene), ***start*** = start position of feature, ***strand*** = forward (+) or reverse (-) strand and more[<sup>6</sup>](#fn6). 

The GTF file is provided in the directory `/mnt/workspace_stud/allstud/homo_sapiens.104.mainChr.gtf`. Some column names are changed to match requirements of `Cicero`. 


```R
gene_anno <- rtracklayer::readGFF("/mnt/workspace_stud/allstud/homo_sapiens.104.mainChr.gtf")

# Change some colnames to match requirements. 
gene_anno$chromosome <- gene_anno$seqid
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name
head(gene_anno)
```


<table class="dataframe">
<caption>A data.frame: 6 × 30</caption>
<thead>
	<tr><th></th><th scope=col>seqid</th><th scope=col>source</th><th scope=col>type</th><th scope=col>start</th><th scope=col>end</th><th scope=col>score</th><th scope=col>strand</th><th scope=col>phase</th><th scope=col>gene_id</th><th scope=col>gene_version</th><th scope=col>⋯</th><th scope=col>exon_number</th><th scope=col>exon_id</th><th scope=col>exon_version</th><th scope=col>protein_id</th><th scope=col>protein_version</th><th scope=col>ccds_id</th><th scope=col>chromosome</th><th scope=col>gene</th><th scope=col>transcript</th><th scope=col>symbol</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>chr1</td><td>havana</td><td>gene      </td><td>11869</td><td>14409</td><td>NA</td><td>+</td><td>NA</td><td>ENSG00000223972</td><td>5</td><td>⋯</td><td>NA</td><td>NA             </td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>chr1</td><td>ENSG00000223972</td><td>NA             </td><td>DDX11L1</td></tr>
	<tr><th scope=row>2</th><td>chr1</td><td>havana</td><td>transcript</td><td>11869</td><td>14409</td><td>NA</td><td>+</td><td>NA</td><td>ENSG00000223972</td><td>5</td><td>⋯</td><td>NA</td><td>NA             </td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>chr1</td><td>ENSG00000223972</td><td>ENST00000456328</td><td>DDX11L1</td></tr>
	<tr><th scope=row>3</th><td>chr1</td><td>havana</td><td>exon      </td><td>11869</td><td>12227</td><td>NA</td><td>+</td><td>NA</td><td>ENSG00000223972</td><td>5</td><td>⋯</td><td>1 </td><td>ENSE00002234944</td><td>1 </td><td>NA</td><td>NA</td><td>NA</td><td>chr1</td><td>ENSG00000223972</td><td>ENST00000456328</td><td>DDX11L1</td></tr>
	<tr><th scope=row>4</th><td>chr1</td><td>havana</td><td>exon      </td><td>12010</td><td>12057</td><td>NA</td><td>+</td><td>NA</td><td>ENSG00000223972</td><td>5</td><td>⋯</td><td>1 </td><td>ENSE00001948541</td><td>1 </td><td>NA</td><td>NA</td><td>NA</td><td>chr1</td><td>ENSG00000223972</td><td>ENST00000450305</td><td>DDX11L1</td></tr>
	<tr><th scope=row>5</th><td>chr1</td><td>havana</td><td>transcript</td><td>12010</td><td>13670</td><td>NA</td><td>+</td><td>NA</td><td>ENSG00000223972</td><td>5</td><td>⋯</td><td>NA</td><td>NA             </td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>chr1</td><td>ENSG00000223972</td><td>ENST00000450305</td><td>DDX11L1</td></tr>
	<tr><th scope=row>6</th><td>chr1</td><td>havana</td><td>exon      </td><td>12179</td><td>12227</td><td>NA</td><td>+</td><td>NA</td><td>ENSG00000223972</td><td>5</td><td>⋯</td><td>2 </td><td>ENSE00001671638</td><td>2 </td><td>NA</td><td>NA</td><td>NA</td><td>chr1</td><td>ENSG00000223972</td><td>ENST00000450305</td><td>DDX11L1</td></tr>
</tbody>
</table>



The gene activity scores are used to get a sense of the accessibility of a promoter and its associated distal elements. 
The function to generate these scores is called `build_gene_activity_matrix`. Its input is a CDS and a Cicero connection list. The output then is an unnormalized table of activity scores of genes. 

In order to create this score a column in the fData table of the input CDS must be added. This column indicates if a peak is a promoter or if the peak is distal `( = "NA")`. 


```R
# Create a gene annotation set that only marks the transcription start sites of the genes. 

# We use this as a proxy for promoters. To do this we need the first exon of each transcript
pos <- subset(gene_anno, strand == "+")
pos <- pos[order(pos$start),] 

# remove all but the first exons per transcript
pos <- pos[!duplicated(pos$transcript),]

# make a 1 base pair marker of the TSS
pos$end <- pos$start + 1 

# Repeat this process for the - strand
neg <- subset(gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 

# remove all but the first exons per transcript
neg <- neg[!duplicated(neg$transcript),] 
neg$start <- neg$end - 1

gene_annotation_sub <- rbind(pos, neg)

# Make a subset of the TSS annotation columns containing just the coordinates and the gene name
gene_annotation_sub <- gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]

# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] <- "gene"

# Add the column with annotated sites to the cds_object
cds_object <- annotate_cds_by_site(cds_object, gene_annotation_sub)
```

Extract the information of the annotated sites of the ***cds_object*** to the variable ***sites***. After changing the format to a data frame, the unused columns can be discarded. Only the information ***gene*** and ***site_name*** (thus peak) is important for the purposes of this work package. 

In the column ***gene*** promoters are called by name and distal elements are named `NA`. We filter this data frame by promoter, so every row with `gene = NA` can be discarded. In the final data frame ***sites*** only promoters and their associated peaks are displayed. 


```R
# Extract data
sites <- fData(cds_object)
sites <- data.frame(sites)

# discard unused columns.
sites <- sites[, c("site_name", "gene")]

# only display peaks that have associated promoters. 
site <- subset(sites, gene != "NA")
rownames(site) <- NULL
head(site)
```


<table class="dataframe">
<caption>A data.frame: 6 × 2</caption>
<thead>
	<tr><th></th><th scope=col>site_name</th><th scope=col>gene</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>chr1_817141_817541</td><td>FAM87B </td></tr>
	<tr><th scope=row>2</th><td>chr1_925472_925872</td><td>SAMD11 </td></tr>
	<tr><th scope=row>3</th><td>chr1_940684_941084</td><td>SAMD11 </td></tr>
	<tr><th scope=row>4</th><td>chr1_970716_971116</td><td>PLEKHN1</td></tr>
	<tr><th scope=row>5</th><td>chr1_973139_973539</td><td>PLEKHN1</td></tr>
	<tr><th scope=row>6</th><td>chr1_982026_982426</td><td>PERM1  </td></tr>
</tbody>
</table>



This data frame is then stored as .csv file in the directory `/mnt/workspace_stud/stud10/cds_sites.csv`.


```R
write.csv(site, '/mnt/workspace_stud/stud10/output/esophagus_mucosa/cds_sites.csv')
```

For demonstration purposes the data is pre-filtered. Instead of calculating the scores for all chromosomes afterwards in `run_cicero()`, we only take those chromosomes that have a promoter. In order to do so a data frame was created to display all chromosomes where promoters were found. 

The variable ***chroms*** is a list of the 8 most  most highly represented chromosomes when counting promoters per chromosome.


```R
# split column to find to find the most represented chromosomes
site[c("chr", "bp1", "bp2")] <- str_split_fixed(site$site_name, '_', 3)

# f is a variable storing occurences of the column chr and order f ascending
f <- table(site$chr)
f <- f[order(f)]

# only use the 8 most highly represented chromosomes, store this list in chroms
f <- as.data.frame(f)
f <- tail(f, n = 8)
chroms <- as.character(f$Var1)
chroms
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'chr7'</li><li>'chr12'</li><li>'chr19'</li><li>'chr3'</li><li>'chr17'</li><li>'chr11'</li><li>'chr2'</li><li>'chr1'</li></ol>



The easiest way to get Cicero co-accessibility scores is to run `run_cicero()`. To do so, you need a Cicero CDS object (created above) and a genome coordinates file, which contains the lengths of each of the chromosomes in your organism. This data is obtained from the .fai file provided in the directory `/mnt/workspace_stud/allstud/`.

The FAI File is an index to enable random access to FASTA and FASTQ Files[<sup>7</sup>](#fn7). Its a text file consisting of TAB-delimited columns for a FASTA file. 
The only information needed for this workpackage is the total length of the reference sequence in bases, an information to be found in column V1 ***chromosome*** and V2 ***length of this sequence***.

To get a total picture of the co-accessibility-scores of the whole genome, all chromosomes must be taken into account. This run is performed with the 8 Chromosomes, obtained from the previous generated variable ***chroms*** and stored in the variable ***chr_len*** displaying these chromosomes with corresponding lengths. 


```R
chr_len <- read.table("/mnt/workspace_stud/allstud/homo_sapiens.104.mainChr.fa.fai", sep = '\t')

# create table with chromosomes matching `chorms`. 
chr_len <- chr_len[chr_len$V1 %in% chroms, c("V1", "V2")]
chr_len
```


<table class="dataframe">
<caption>A data.frame: 8 × 2</caption>
<thead>
	<tr><th></th><th scope=col>V1</th><th scope=col>V2</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>chr1 </td><td>248956422</td></tr>
	<tr><th scope=row>2</th><td>chr2 </td><td>242193529</td></tr>
	<tr><th scope=row>3</th><td>chr3 </td><td>198295559</td></tr>
	<tr><th scope=row>7</th><td>chr7 </td><td>159345973</td></tr>
	<tr><th scope=row>11</th><td>chr11</td><td>135086622</td></tr>
	<tr><th scope=row>12</th><td>chr12</td><td>133275309</td></tr>
	<tr><th scope=row>17</th><td>chr17</td><td> 83257441</td></tr>
	<tr><th scope=row>19</th><td>chr19</td><td> 58617616</td></tr>
</tbody>
</table>



The wrapper function `run_cicero()` runs the primary functions of the Cicero pipeline. It takes the previously generated CDS object using `make_cicero_cds()` as input. Furthermore this function needs genomic coordinations, a data frame containing chromosome lengths in base pairs. The argument silent is whether to print progress messages and the sample_num is how many sample genomic windows to use to generate the distance parameter. The distance parameter can be called in a separated function (`estimate_distance_parameter()`). 


```R
conns <- run_cicero(cicero_cds, 
                    chr_len,  
                    sample_num = 100)
head(conns)
```

    [1] "Starting Cicero"
    [1] "Calculating distance_parameter value"
    [1] "Running models"
    [1] "Assembling connections"
    [1] "Successful cicero models:  4812"
    [1] "Other models: "
    
    Zero or one element in range 
                             229 
    [1] "Models with errors:  0"
    [1] "Done"



<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>Peak1</th><th scope=col>Peak2</th><th scope=col>coaccess</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>chr11_100004558_100004958</td><td>chr11_99929877_99930277  </td><td>0</td></tr>
	<tr><th scope=row>2</th><td>chr11_100004558_100004958</td><td>chr11_99988198_99988598  </td><td>0</td></tr>
	<tr><th scope=row>3</th><td>chr11_100004558_100004958</td><td>chr11_99988668_99989068  </td><td>0</td></tr>
	<tr><th scope=row>5</th><td>chr11_100004558_100004958</td><td>chr11_100066418_100066818</td><td>0</td></tr>
	<tr><th scope=row>6</th><td>chr11_100004558_100004958</td><td>chr11_100067321_100067721</td><td>0</td></tr>
	<tr><th scope=row>7</th><td>chr11_100004558_100004958</td><td>chr11_100068804_100069204</td><td>0</td></tr>
</tbody>
</table>



For further usage the table of co-accessibility-scores is stored as .csv in the directory `/mnt/workspace_stud/stud10/`.


```R
write.csv(conns, '/mnt/workspace_stud/stud10/output/esophagus_mucosa/conns.csv')
```

`Cicero`holds a function to calculate the so called "gene activity score". 
The resulting matrix created by `build_gene_activity_matrix()` is unnormalized and must be normalized by the function `normalize_gene_activites()`. To call this function a vector of total accessible sites per cell has to be provided. This vector can be found in the pData table of the CDS ***num_genes_expressed***. 


```R
#### Generate gene activity scores ####
# generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(cds_object, conns)

# remove any rows/columns with all zeroes
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                       !Matrix::colSums(unnorm_ga) == 0]

# make a list of num_genes_expressed
num_genes <- pData(cds_object)$num_genes_expressed
names(num_genes) <- row.names(pData(cds_object))

# normalize
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)

# if you had two datasets to normalize, you would pass both:
# num_genes should then include all cells from both sets
unnorm_ga2 <- unnorm_ga
cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga, unnorm_ga2), 
                                                    num_genes)
```

The activity score matrix is stored in a list, and must be transformed to a matrix by `unlist()`. 


```R
cicero_gene_activities <- unlist(cicero_gene_activities[[1]])
```

For further usage the table of gene activities is stored as .mtx in the directory `/mnt/workspace_stud/stud10/`.
Later on the Matrix will be assigned with the clusters as column names, not the barcodes. So only the rownames has to be stored with the ***cicero_gene_activities***-Matrix. 


```R
write(x = rownames(cicero_gene_activities), 
      file = "/mnt/workspace_stud/stud10/output/esophagus_mucosa/cicero_gene_activities_rows.tsv")
write(x = colnames(cicero_gene_activities), 
      file = "/mnt/workspace_stud/stud10/output/esophagus_mucosa/cicero_gene_activities_cols.tsv")
writeMM(cicero_gene_activities, '/mnt/workspace_stud/stud10/output/esophagus_mucosa/cicero_genes_activity.mtx')
```


    NULL


***
# Appendix

### References 

<sup>1</sup><span id="fn1"> **Introduction to ATAC-seq:** </span>  Sun, Yuanyuan, Nan Miao, and Tao Sun. "Detect accessible chromatin using ATAC-sequencing, from principle to applications." Hereditas 156.1 (2019): 1-9.

<sup>2</sup><span id="fn2"> **Origin study:**</span> Zhang, Kai, et al. "A cell atlas of chromatin accessibility across 25 adult human tissues." BioRxiv (2021).

<sup>3</sup><span id="fn3"> **Cicero:**</span> Pliner, Hannah A., et al. "Cicero predicts cis-regulatory DNA interactions from single-cell chromatin accessibility data." Molecular cell 71.5 (2018): 858-871.

<sup>4</sup><span id="fn4"> **Cicero Code:**</span>  https://github.com/cole-trapnell-lab/cicero-release/tree/master/man

<sup>5</sup><span id="fn5"> **UMAP:**</span> https://pair-code.github.io/understanding-umap/

<sup>6</sup><span id="fn6"> **GTF File Format**</span>  http://www.ensembl.org/info/website/upload/gff.html

<sup>7</sup><span id="fn7"> **FAI File Format**</span> **FAI File Format:** http://www.htslib.org/doc/faidx.html

****
### List of Figures 

- Figure 1: https://www.cell.com/molecular-cell/fulltext/S1097-2765(18)30547-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1097276518305471%3Fshowall%3Dtrue#figures
