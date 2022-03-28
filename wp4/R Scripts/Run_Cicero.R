library(ggplot2)
library(dplyr)
library(monocle3)
library(cicero)
library(seqTools)
library(stringr)
library(Matrix)

gene_anno <- rtracklayer::readGFF("/mnt/workspace_stud/allstud/homo_sapiens.104.mainChr.gtf")

# Change some colnames to match requirements. 
gene_anno$chromosome <- gene_anno$seqid
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name
head(gene_anno)

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

# Extract data
sites <- fData(cds_object)
sites <- data.frame(sites)

# discard unused columns.
sites <- sites[, c("site_name", "gene")]

# only display peaks that have associated promoters. 
site <- subset(sites, gene != "NA")
rownames(site) <- NULL
head(site)

write.csv(site, '/mnt/workspace_stud/stud10/output/esophagus_mucosa/cds_sites.csv')

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

chr_len <- read.table("/mnt/workspace_stud/allstud/homo_sapiens.104.mainChr.fa.fai", sep = '\t')

# create table with chromosomes matching `chorms`. 
chr_len <- chr_len[chr_len$V1 %in% chroms, c("V1", "V2")]
chr_len

conns <- run_cicero(cicero_cds, 
                    chr_len,  
                    sample_num = 100)
head(conns)

write.csv(conns, '/mnt/workspace_stud/stud10/output/esophagus_mucosa/conns.csv')

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

cicero_gene_activities <- unlist(cicero_gene_activities[[1]])
cicero_gene_activities

write(x = rownames(cicero_gene_activities), 
      file = "/mnt/workspace_stud/stud10/output/esophagus_mucosa/cicero_gene_activities_rows.tsv")
write(x = colnames(cicero_gene_activities), 
      file = "/mnt/workspace_stud/stud10/output/esophagus_mucosa/cicero_gene_activities_cols.tsv")
writeMM(cicero_gene_activities, '/mnt/workspace_stud/stud10/output/esophagus_mucosa/cicero_genes_activity.mtx')