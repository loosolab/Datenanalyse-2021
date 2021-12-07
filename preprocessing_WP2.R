library(SnapATAC);
y.sp = createSnap(
  file="/home/rstudio/workspaces/stud4/SnaptoolsTest/ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver.snap",
  sample="ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver",
  num.cores=5
);
print(y.sp)

showBinSizes("/home/rstudio/workspaces/stud4/SnaptoolsTest/ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver.snap");
y.sp = addBmatToSnap(y.sp, bin.size=5000, num.cores=5);
y.sp = makeBinary(y.sp, mat="bmat");

y.sp

bin.cov = log10(Matrix::colSums(y.sp@bmat)+1);
bin.cov

hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 5)
);

bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
y.sp = y.sp[, idy, mat="bmat"];
y.sp


y.sp = runDiffusionMaps(
  obj=y.sp,
  input.mat="bmat", 
  num.eigs=50
);

plotDimReductPW(
  obj=y.sp, 
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

y.sp = runKNN(
  obj=y.sp,
  eigs.dims=1:7,
  k=15
);

y.sp=runCluster(
  obj=y.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  seed.use=10
);

y.sp@metaData$cluster = y.sp@cluster;

y.sp = runViz(
  obj=y.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:7, 
  method="Rtsne",
  seed.use=10
);

par(mfrow = c(2, 2));

plotViz(
  obj=y.sp,
  method="tsne", 
  main="Right lobe of liver Cluster",
  point.color=y.sp@cluster, 
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

y.sp@barcode
ca_table = do.call(rbind, Map(data.frame, cell_name=y.sp@barcode, cluster=y.sp@cluster))
write.table(ca_table,"cluster_assignment_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

ca_table