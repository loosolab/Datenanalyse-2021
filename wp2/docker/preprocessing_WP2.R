library(GenomicRanges)
library(SnapATAC)
library(viridisLite)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
print("base_path ist:")
print(args[1])
print("name is :")
print(args[2])

genefile = paste(args[1],"gencode.filtered.bed",sep="")
print(genefile)
blacklist_regions = paste(args[1],"hg38-blacklist.v2_parsed.bed",sep="")
sample_name = "args[2]"

setwd(args[1])

# load .snap file
s_file=paste(args[1],args[2],"/",args[2],".snap",sep="")
#s_file="/home/rstudio/workspaces/stud4/SnaptoolsTest/ENC-1JKYN-020-SM-C1PX3_snATAC_thoracic_aorta.snap"
x.sp = createSnap(
  file=s_file,
  sample=sample_name,
  num.cores=6
);
print("snap createt")
# add bin matrix
x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=6);
x.sp = makeBinary(x.sp, mat="bmat");
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
print("chroms removed")
# calc log10(bin coverage)
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
jpeg(file=paste(args[2],"_coverage_hist.jpeg",sep=""))
# plot histogram of bin coverage
hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 5)
);
dev.off()

#system('touch /mnt/testing/saving_plot1.jpg')
#jpeg(file="saving_plot1.jpeg")
#dev.off()
print("hist done")
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
print("save file in:")
print(args[1])
print("name:")
print(paste(args[2],"_Eigenvector_plot.jpeg"),sep="")


jpeg(file=paste(args[2],"_Eigenvector_plot.jpeg",sep=""))

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

print("vectors done")
# select first pair that looks like a blob
#v_pair = 5

# continue dimensional reduction with first n eigenvectors
#x.sp = runKNN(
#  obj=x.sp,
#  eigs.dims=1:v_pair, # choose first plot with "blob"
#  k=15
#);

# perform clustering
#x.sp=runCluster(
#  obj=x.sp,
#  tmp.folder=tempdir(),
#  louvain.lib="R-igraph",
#  seed.use=10
#);

#x.sp@metaData$cluster = x.sp@cluster;

# calculate t-sne
#x.sp = runViz(
#  obj=x.sp, 
#  tmp.folder=tempdir(),
#  dims=2,
#  eigs.dims=1:v_pair, 
#  method="Rtsne",
#  seed.use=10
#);

# plot t-sne
#plotViz(
#  obj=x.sp,
#  method="tsne", 
#  main=sample_name,
#  point.color=x.sp@cluster, 
#  point.size=1, 
#  point.shape=19, 
#  point.alpha=0.8, 
#  text.add=TRUE,
#  text.size=1.5,
#  text.color="black",
#  text.halo.add=TRUE,
#  text.halo.color="white",
#  text.halo.width=0.2,
#  down.sample=10000,
#  legend.add=FALSE
#);

## new gene annotation
## TODO
## # gene annotation
## genes = read.table(genefile);
## genes.gr = GRanges(genes[,1], 
##                      IRanges(genes[,2], genes[,3]), name=genes[,4]
## );
## 
## x.sp = addBmatToSnap(x.sp);
## # TODO: Select highest "expressed" genes of each cluster
## x.sp = createGmatFromMat(
##   obj=x.sp, 
##   input.mat="bmat",
##   genes=genes.gr,
##   do.par=TRUE,
##   num.cores=6
## );

## generate cluster assignment table

#setwd("/mnt/testing")
#ca_table = do.call(rbind, Map(data.frame, cell_name=x.sp@barcode, cluster=x.sp@cluster))
#write.table(ca_table,"cluster_assignment_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#print("cluster table done")

# call peaks for all cluster with more than 200 cells
#clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 200)];
#peaks.ls = mclapply(seq(clusters.sel), function(i){
#  print(clusters.sel[i]);
#  runMACS(
#    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),],
#    output.prefix=paste0(paste(sample_name, ".", sep=""), gsub(" ", "_", clusters.sel)[i]),
#    path.to.snaptools="/usr/local/bin/snaptools",
#    path.to.macs="/usr/local/bin/macs2",
#    gsize="hs", # mm, hs, etc
#    buffer.size=500, 
#    num.cores=6,
#    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
#    tmp.folder=tempdir()
#  );
#}, mc.cores=6);
# assuming all .narrowPeak files in the current folder are generated from the clusters
#peaks.names = system("ls | grep narrowPeak", intern=TRUE);
#peak.gr.ls = lapply(peaks.names, function(x){
#  peak.df = read.table(x)
#  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
#})
#peak.gr = reduce(Reduce(c, peak.gr.ls));
#peak.gr

# create a cell-by-peak matrix
#peaks.df = as.data.frame(peak.gr)[,1:3];
#write.table(peaks.df,file = "peaks.combined.bed",append=FALSE,
#              quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
#              row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
#              fileEncoding = "")

#saveRDS(x.sp, file=paste(sample_name, ".rds", sep=""))

# TODO: don't think that the commented steps are necessary for our pipeline?
# create cell-by-peak matrix and add to the snap file
# Terminal: snaptools snap-add-pmat \
# --snap-file /home/rstudio/workspaces/stud4/SnaptoolsTest/ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver.snap \
# --peak-file peaks.combined.bed

# add cell-by-peak matrix
# x.sp = readRDS("right_lobe_of_liver.snap.rds");
# x.sp = addPmatToSnap(x.sp);
# x.sp = makeBinary(x.sp, mat="pmat");

#?paste

#ptest = paste(sample_name, ".", sep="")
#ptest

#x.sp
