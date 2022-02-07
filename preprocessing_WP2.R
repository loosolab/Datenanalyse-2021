library(GenomicRanges)
library(SnapATAC)
library(viridisLite)
library(ggplot2)
library(rtracklayer)

genefile = "/home/rstudio/workspaces/stud4/genefile/gencode.filtered.bed"
blacklist_regions = "/home/rstudio/workspaces/stud4/SnaptoolsTest/blacklists/hg38-blacklist.v2_parsed.bed"

# load .snap file
s_file="/home/rstudio/workspaces/stud4/SnaptoolsTest/ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver.snap"
sample_name = "right-lobe-of-liver"
s_file="/home/rstudio/workspaces/stud4/SnaptoolsTest/ENC-1JKYN-020-SM-C1PX3_snATAC_thoracic_aorta.snap"
sample_name = "thoracic-aorta"

# set/getwd and create folder
setwd("/home/rstudio/workspaces/stud3/rstudio_workspace/Datenanalyse_2/Datenanalyse-2021")
dir.create(file.path(getwd(), sample_name))
setwd(file.path(getwd(), sample_name))

x.sp = createSnap(
  file=s_file,
  sample=sample_name,
  num.cores=6
);
x.sp

# add bin matrix
x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=6);
x.sp = makeBinary(x.sp, mat="bmat");
x.sp

# barcode filtering
gtf.gr <- rtracklayer::import("/home/rstudio/workspaces/mbentse/homo_sapiens/homo_sapiens.104.mainChr.gtf")
gene.gr <- gtf.gr[gtf.gr$type == "gene"]
# extract promoter region for each gene
promoter.gr <- reduce(promoters(gene.gr, upstream=2000, downstream=0))
# find promoter overlapping bins 
ov = findOverlaps(x.sp@feature, promoter.gr);
idy = queryHits(ov);
log_cov = log10(SnapATAC::rowSums(x.sp, mat="bmat")+1)
promoter_ratio = Matrix::rowSums(x.sp@bmat[,idy]) / Matrix::rowSums(x.sp@bmat);
png(filename="barcodes.png", width=6, height=4, units="in", res=300)
plot(log_cov, promoter_ratio, cex=0.5, col="grey", xlab="log(count(bins))", ylab="Promoter ratio", ylim=c(0,1 ));
dev.off()

# chose the cutoff based on the plot 
# liver 0.3 0.6 2.5
idx = which(promoter_ratio > 0.3 & promoter_ratio < 0.6 & log_cov > 2.5);

x.sp = x.sp[idx,]
x.sp

# bin filtering
black_list = read.table(blacklist_regions);
black_list.gr = GRanges(
  black_list[,1], 
  IRanges(black_list[,2], black_list[,3])
);
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr));
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp

# remove unwanted chromosomes
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp

# calc log10(bin coverage)
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);

# plot histogram of bin coverage
png(filename="hist.png", width=6, height=4, units="in", res=300)
hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 5)
);
dev.off()


# filter bins with low coverage
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];
x.sp

# dimensional reduction
x.sp = runDiffusionMaps(
  obj=x.sp,
  input.mat="bmat", 
  num.eigs=50
);

# show first 25 pairs of eigenvectors
png(filename="eigenvectors.png", width=6, height=4, units="in", res=300)
plotDimReductPW(
  obj=x.sp, 
  eigs.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=5000,
  pdf.file.name=NULL, 
  pdf.height=7, 
  pdf.width=7
);
dev.off()
# select first pair that looks like a blob
# liver: 7
# aorta:
v_pair = 7

# continue dimensional reduction with first n eigenvectors
x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:v_pair, # choose first plot with "blob"
  k=15
);

# perform clustering
x.sp=runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  seed.use=10
);

x.sp@metaData$cluster = x.sp@cluster;

# calculate t-sne
x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:v_pair, 
  method="Rtsne",
  seed.use=10
);

# plot t-sne
png(filename="t-sne.png", width=6, height=6, units="in", res=300)
plotViz(
  obj=x.sp,
  method="tsne", 
  main=sample_name,
  point.color=x.sp@cluster, 
  point.size=0.8, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  #down.sample=10000,
  legend.add=FALSE
);
dev.off()

# new gene annotation
# TODO
# # gene annotation
# genes = read.table(genefile);
# genes.gr = GRanges(genes[,1], 
#                      IRanges(genes[,2], genes[,3]), name=genes[,4]
# );
# 
# x.sp = addBmatToSnap(x.sp);
# # TODO: Select highest "expressed" genes of each cluster
# x.sp = createGmatFromMat(
#   obj=x.sp, 
#   input.mat="bmat",
#   genes=genes.gr,
#   do.par=TRUE,
#   num.cores=6
# );

# generate cluster assignment table
ca_table = do.call(rbind, Map(data.frame, cell_name=x.sp@barcode, cluster=x.sp@cluster))
write.table(ca_table,"cluster_assignment_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# call peaks for all cluster with more than 100 cells
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 100)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0(paste(sample_name, ".", sep=""), gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/usr/local/bin/snaptools",
    path.to.macs="/usr/local/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=6,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir()
  );
}, mc.cores=6);
# assuming all .narrowPeak files in the current folder are generated from the clusters
peaks.names = system("ls | grep narrowPeak", intern=TRUE);
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls));
peak.gr

# create a cell-by-peak matrix
peaks.df = as.data.frame(peak.gr)[,1:3];
write.table(peaks.df,file = "peaks.combined.bed",append=FALSE,
              quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")

saveRDS(x.sp, file=paste(sample_name, ".rds", sep=""))

# TODO: don't think that the commented steps are necessary for our pipeline?
# create cell-by-peak matrix and add to the snap file
# Terminal: snaptools snap-add-pmat \
# --snap-file /home/rstudio/workspaces/stud4/SnaptoolsTest/ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver.snap \
# --peak-file peaks.combined.bed

# add cell-by-peak matrix
# x.sp = readRDS("right_lobe_of_liver.snap.rds");
# x.sp = addPmatToSnap(x.sp);
# x.sp = makeBinary(x.sp, mat="pmat");