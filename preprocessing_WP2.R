library(GenomicRanges)
library(SnapATAC);
genefile = "/home/rstudio/workspaces/stud4/genefile/gencode.hg38.gene.bed"
# use  one .snap file only...
s_file="/home/rstudio/workspaces/stud4/SnaptoolsTest/ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver.snap"
x.sp = createSnap(
  file=s_file,
  sample="ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver",
  num.cores=6
);
x.sp

# ...or merge .snap files first
file.list = c("/home/rstudio/workspaces/stud4/SnaptoolsTest/ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver.snap", 
              "/home/rstudio/workspaces/stud4/SnaptoolsTest/ENC-1JKYN-020-SM-C1PX3_snATAC_thoracic_aorta.snap");
sample.list = c("liver", "aorta");
x.sp = createSnap(file=file.list, sample=sample.list);
x.sp

# add bin matrix
# showBinSizes(s_file);
x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=6);
x.sp = makeBinary(x.sp, mat="bmat");
x.sp

# calc log10(bin coverage)
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);

# plot histogram of bin coverage
hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 5)
);

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

# select first pair that looks like a blob
v_pair = 10

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

par(mfrow = c(2, 2));

# plot t-sne
plotViz(
  obj=x.sp,
  method="tsne", 
  main="Right lobe of liver",
  point.color=x.sp@cluster, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
);

# gene annotation
genes = read.table(genefile);
genes.gr = GRanges(genes[,1], 
                     IRanges(genes[,2], genes[,3]), name=genes[,4]
);
marker.genes = c(
  "Snap25", "Gad2", "Apoe",
  "C1qb", "Pvalb", "Vip", 
  "Sst", "Lamp5", "Slc17a7"
);
genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)];
# re-add the cell-by-bin matrix to the snap object;
x.sp = addBmatToSnap(x.sp);
x.sp = createGmatFromMat(
  obj=x.sp, 
  input.mat="bmat",
  genes=genes.sel.gr,
  do.par=TRUE,
  num.cores=10
);
# normalize the cell-by-gene matrix
x.sp = scaleCountMatrix(
  obj=x.sp, 
  cov=x.sp@metaData$passed_filters + 1,
  mat="gmat",
  method = "RPM"
);
# smooth the cell-by-gene matrix
x.sp = runMagic(
  obj=x.sp,
  input.mat="gmat",
  step.size=3
);
par(mfrow = c(3, 3));
for(i in 1:9){
  plotFeatureSingle(
    obj=x.sp,
    feature.value=x.sp@gmat[, marker.genes[i]],
    method="tsne", 
    main=marker.genes[i],
    point.size=0.1, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0, 1)
  )};

# generate cluster assignment table
x.sp@barcode
ca_table = do.call(rbind, Map(data.frame, cell_name=x.sp@barcode, cluster=x.sp@cluster))
write.table(ca_table,"cluster_assignment_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

ca_table