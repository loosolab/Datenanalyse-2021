library(Matrix)
library(stringr)
library(monocle3)
library(cicero)

indata <- Matrix::readMM("/mnt/workspace_stud/stud2/output/esophagus_mucosa/wp4/esophagus_mucosa_X.mtx")
indata@x[indata@x > 0] <- 1
indata <- t(indata)
# filter out peaks with less than 30 cells 
indata <- indata[rowSums(indata) > 30,]
head(indata)

cellinfo <- read.csv("/mnt/workspace_stud/stud2/output/esophagus_mucosa/wp4/esophagus_mucosa_obs.csv", 
                     sep = '+', row.names = NULL)
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

bar_cluster <- paste(cellinfo$barcode, cellinfo$cluster, sep="_")
write.csv(bar_cluster, '/mnt/workspace_stud/stud10/output/esophagus_mucosa/bar_cluster.tsv')

# format peak info
peakinfo <- read.csv("/mnt/workspace_stud/stud2/output/esophagus_mucosa/wp4/esophagus_mucosa_var.csv",
                     header = TRUE, sep = '_', row.names = NULL, col.names= c("chr", "bp1", "shared", "cells.variability", "score"))

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

# Creating new column and paste this concatenated string as new row names. 
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")

# Set rownames of the dataframe
row.names(peakinfo) <- peakinfo$site_name

# Discard unused columns. 
peakinfo <- peakinfo[, c("site_name", "chr", "bp1", "bp2")]
head(peakinfo)

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

cds_object <- suppressWarnings(new_cell_data_set(indata,
                                                 cell_metadata = cellinfo, 
                                                 gene_metadata = peakinfo))

cds_object <- monocle3::detect_genes(cds_object)
cds_object <- cds_object[Matrix::rowSums(exprs(cds_object)) != 0,] 
