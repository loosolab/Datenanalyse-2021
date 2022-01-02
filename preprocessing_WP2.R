library(SnapATAC);
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

# plot t-sne
x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:v_pair, 
  method="Rtsne",
  seed.use=10
);

par(mfrow = c(2, 2));

# generate cluster assignment table
x.sp@barcode
ca_table = do.call(rbind, Map(data.frame, cell_name=x.sp@barcode, cluster=x.sp@cluster))
write.table(ca_table,"cluster_assignment_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

ca_table