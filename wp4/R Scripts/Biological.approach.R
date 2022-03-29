library(dplyr)
library(stringr)
library(Matrix)


######## Research Question 1 ########
conns <- read.csv("/mnt/workspace_stud/stud10/output/esophagus_mucosa/conns.csv", 
                  sep = ',', row.names = NULL)
conns <- conns[c("Peak1", "Peak2", "coaccess")]

cds_site <- read.csv("/mnt/workspace_stud/stud10/output/esophagus_mucosa/cds_sites.csv", 
                     sep = ',', row.names = NULL)
colnames(cds_site) <- c("num", "Peak2", "gene")
cds_site <- cds_site[,c("Peak2", "gene")]

# combine co-access table with cds site table to assign promoters to peaks
dff <- inner_join(conns, cds_site, by = "Peak2")
# remove entries with NA
dff <- dff[complete.cases(dff),]
# order data by gene with inner order by coaccess
dff <- dff[order(dff$gene, dff$coaccess),]
dff <- subset(dff, coaccess > 0.2 | coaccess < -0.01)

oc_dff <- dff
oc_dff[oc_dff < 0] <- '-'
oc_dff$coaccess[oc_dff$coaccess > 0] <- '+'

list_genes <-sort(table(unlist(oc_dff$gene)))
list_genes <- as.data.frame(list_genes)
colnames(list_genes) <- c("gene", "Freq")
list_genes <- list_genes[list_genes$Freq>30, ]
gene_dff <- oc_dff[which(oc_dff$gene == 'ST6GAL1'), ]

cluster_peaks <- read.csv("/mnt/workspace_stud/stud10/cluster_peaks.csv", 
                          sep = ',', row.names = NULL)
cluster_peaks <- cluster_peaks[complete.cases(cluster_peaks), ]
colnames(cluster_peaks) <- c("Peak2", "0", "1", "2", "3", "4", "5", 
                             "6", "7", "8", "9", "10", "11", "12",
                             "13", "14", "15", "16", "17")

cluster_peaks$max <- c("9_2", "9", "13_20", "2", "13")
cluster_peaks <- cluster_peaks[, c("Peak2", "max")]

cluster_gene_dff <- inner_join(gene_dff, cluster_peaks, by = "Peak2")
cluster_gene_dff <- cluster_gene_dff[order(cluster_gene_dff$Peak2, 
                                           cluster_gene_dff$coaccess),]

list_max_values <- unique(cluster_gene_dff$max)
len <- length(list_max_values)

lhs  <- paste("sub",    1:len,     sep="")
rhs  <- paste("list_max_values[[",1:len,"]]", sep="")
eq   <- paste(paste(lhs, rhs, sep="<-"), collapse=";")
eval(parse(text=eq))

lhs  <- paste("sublist",    1:len,     sep="")
rhs  <- paste("cluster_gene_dff[which(cluster_gene_dff$max == sub",1:len, "),]", sep="")
eq   <- paste(paste(lhs, rhs, sep="<-"), collapse=";")
eval(parse(text=eq))
# create list with sublists 
sublist <- list(sublist1, sublist2, sublist3, sublist4, sublist5)

counter1 <- 0
counter2 <- 0
for (a in sublist){
  counter1 <<- counter1 + 1
  for (b in sublist){
    counter2 <<- counter2 + 1
    ij_list <- inner_join(a, b, by = "Peak1")
    ij_list <- ij_list[, c("Peak1", "coaccess.x", "max.x", "coaccess.y", "max.y")]
    x <- ij_list$coaccess.x
    y <- ij_list$coaccess.y
    if (! identical(x, y)){
      print(paste(counter1, "not equal to ", counter2))
    }
    else {
    }
    
  }
  counter2 <- 0
  print(paste(counter1, " is equal with every data frame."))
}

freq <- as.data.frame(table(cluster_gene_dff$Peak1))
freq <- freq[which(freq$Freq == 1), ]
colnames(freq) <- c("Peak1", "Freq")

cluster_gene_dff_unique <- inner_join(freq, cluster_gene_dff, by = "Peak1")
cluster_gene_dff_unique <- cluster_gene_dff_unique[order(cluster_gene_dff_unique$max, 
                                                         cluster_gene_dff_unique$coaccess),]

# create data frame with requirements to filter cl_indata 
coacc_dff <- gene_dff
rownames(coacc_dff) <- NULL
coacc_dff <- coacc_dff[!duplicated(coacc_dff$Peak2),]
row.names(coacc_dff) <- coacc_dff$Peak2

# store filtered matrix in directory
write(x = rownames(subset), 
      file = "/mnt/workspace_stud/stud10/output/esophagus_mucosa/peaks.tsv")
write(x = colnames(subset), 
      file = "/mnt/workspace_stud/stud10/output/esophagus_mucosa/cluster.tsv")

subset <- Matrix(subset , sparse = T )
writeMM(obj = subset, 
        file = "/mnt/workspace_stud/stud10/output/esophagus_mucosa/cl_indata_esophagus_mucosa.mtx")

act <- Matrix::readMM("/mnt/workspace_stud/stud10/output/esophagus_mucosa/cicero_genes_activity.mtx")
rows <- read.csv('/mnt/workspace_stud/stud10/output/esophagus_mucosa/cicero_gene_activities_rows.tsv', 
                 col.names = 'rows', sep = ',', row.names = NULL, header = FALSE)
cols <- read.csv('/mnt/workspace_stud/stud10/output/esophagus_mucosa/cicero_gene_activities_cols.tsv',  
                 col.names = 'barcode', sep = ',', header = FALSE)

bar_cluster <- read.csv('/mnt/workspace_stud/stud10/output/esophagus_mucosa/bar_cluster.tsv',  
                        col.names = c("num", 'barcode'), sep = ',', header = FALSE)
bar_cluster[c('barcode', 'cluster')] <- str_split_fixed(bar_cluster$barcode, '_', 2)
bar_cluster_df <- bar_cluster[-1,]
bar_cluster_df$num <- NULL
barcode_dff <- inner_join(cols, bar_cluster_df , by = "barcode")
barcode_dff$joined <- paste(barcode_dff$barcode, barcode_dff$cluster, sep="_")

# set col and rown names of act.
row.names(act) <- rows$rows
colnames(act) <- barcode_dff$joined

include_list <- c("ST6GAL1")
gene_matrix <- suppressWarnings(subset(act, rownames(act) %in% include_list))
gene_matrix <- cbind(gene_matrix, barcode_dff$joined)
colnames(gene_matrix) <- c("activity_score", "barcode_cluster")

gene_matrix <- as.data.frame(gene_matrix)
gene_matrix[c('barcode','cluster')] <- str_split_fixed(gene_matrix$barcode_cluster, '_', 2)
gene_matrix <- gene_matrix[gene_matrix$activity_score > 0.1, ]
gene_matrix <- gene_matrix[, c("activity_score", "cluster")]
gene_matrix <- gene_matrix[order(gene_matrix$activity_score), ]

cluster_max <- gene_matrix[gene_matrix$activity_score == max(gene_matrix$activity_score), ]
cluster <- cluster_max$cluster
cluster

######## Research Question 2 ########
conns_all <- read.csv("/mnt/workspace_stud/stud10/output/conns_all.csv", 
                      sep = ',', row.names = NULL)
conns_all <- conns_all[c("Peak1", "Peak2", "coaccess", "gene", "tissue")]

f <- conns_all %>% 
  count(gene) %>% 
  filter(n == 1) %>% 
  select(-n)

unique_genes <- inner_join(conns_all, f, by = "gene")
unique_genes <- unique_genes[, c("Peak2", "gene", "tissue")]

len <- unique_genes[unique_genes$tissue == 'esophagus_mucosa', ]
length(len$tissue)

######## Cicero Connection Plot ########
library(cicero)
library(ggplot2)
library(stringr)

gene_anno <- rtracklayer::readGFF("/mnt/workspace_stud/allstud/homo_sapiens.104.mainChr.gtf")
# rename some columns to match requirements
gene_anno$chromosome <- gene_anno$seqid
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

filtered_dff <- subset(dff, coaccess > 0.2)
g <- sort(table(filtered_dff$gene))
g <- tail(g, n=1)
g <- names(g)
dff_gene <- filtered_dff[filtered_dff$gene == g,]
dff_gene <- tail(dff_gene, n = 10)
tail(dff_gene)

minbp <- dff_gene[order(dff_gene$Peak1), ]
minbp <- minbp[1, ]
minbp[c('chr', 'bp1', 'bp2')] <- str_split_fixed(minbp$Peak1, '_', 3)
view <- minbp$Peak2
chrom <- minbp$chr
minbp <- as.numeric(minbp$bp1)

maxbp <- dff_gene[order(dff_gene$Peak2), ]
maxbp <- tail(maxbp, n = 1 )
maxbp[c('chr', 'bp1', 'bp2')] <- str_split_fixed(maxbp$Peak2, '_', 3)
maxbp <- as.numeric(maxbp$bp2)

plot_connections(dff, chrom, minbp, maxbp,
                 viewpoint = view,
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.20, 
                 connection_width = 1.5, 
                 collapseTranscripts = "longest" )



