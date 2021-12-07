library(rhdf5);
file = "/home/rstudio/workspaces/stud3/TestData/ENC-1JKYN-012-SM-JF1O6_snATAC_body_of_pancreas2.snap"
idx <- h5read(file, paste("AM", 5000, "idx", sep="/"));
head(idx)

library(SnapATAC);
z.sp = createSnap(
  file="/home/rstudio/workspaces/stud4/SnaptoolsTest/ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver2.snap",
  sample="ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver2",
  num.cores=5
);
print(z.sp)

showBinSizes("/home/rstudio/workspaces/stud4/SnaptoolsTest/ENC-1LGRB-069-SM-A8WNZ_snATAC_right_lobe_of_liver2.snap");
z.sp = addBmatToSnap(z.sp, bin.size=5000, num.cores=5);
z.sp = makeBinary(z.sp, mat="bmat");

z.sp

bin.cov = log10(Matrix::colSums(z.sp@bmat)+1);
bin.cov

hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 5)
);

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


plotFeatureSingle(
  obj=y.sp,
  feature.value=log(y.sp@metaData[,"passed_filters"]+1,10),
  method="tsne", 
  main="Right lobe of liver read depth",
  point.size=0.2, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99)
); 

plotFeatureSingle(
  obj=y.sp,
  feature.value=y.sp@metaData$peak_region_fragments / y.sp@metaData$passed_filters,
  method="tsne", 
  main="Right lobe of liver FRiP",
  point.size=0.2, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99) # remove outliers
);

plotFeatureSingle(
  obj=y.sp,
  feature.value=y.sp@metaData$duplicate / y.sp@metaData$total,
  method="tsne", 
  main="Right lobe of liver Duplicate",
  point.size=0.2, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99) # remove outliers
)
  

library(SnapATAC);
x.sp = createSnap(
  file="/home/rstudio/workspaces/stud3/TestData/ENC-1LVAN-194-SM-JF1NY_snATAC_heart_left_ventricle.demultiplexed.snap",
  sample="artery_aorta_SM-C1PX3",
  num.cores=1
);
print(x.sp)

showBinSizes("/home/rstudio/workspaces/stud3/TestData/ENC-1LVAN-194-SM-JF1NY_snATAC_heart_left_ventricle.demultiplexed.snap");
x.sp = addBmatToSnap(x.sp, bin.size=10000, num.cores=5);
x.sp = makeBinary(x.sp, mat="bmat");

x.sp



bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
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
x.sp = x.sp[, idy, mat="bmat"];
x.sp

x.sp@bmat
x.sp@barcode

x.sp = runDiffusionMaps(
  obj=x.sp,
  input.mat="bmat", 
  num.eigs=50
);


barcodes = read.csv(
  "/home/rstudio/stud3/data/snapATAC_examples/example_1/atac_v1_adult_brain_fresh_5k_singlecell.csv",
  head=TRUE
);
barcodes = barcodes[2:nrow(barcodes),];
promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
UMI = log(barcodes$passed_filters+1, 10);
data = data.frame(UMI=UMI, promoter_ratio=promoter_ratio);
barcodes$promoter_ratio = promoter_ratio;
library(viridisLite);
library(ggplot2);
p1 = ggplot(
  data,
  aes(x= UMI, y= promoter_ratio)) +
  geom_point(size=0.1, col="grey") +
  theme_classic() +
  ggtitle("10X Fresh Adult Brain") +
  ylim(0, 1) + xlim(0, 6) +
  labs(x = "log10(UMI)", y="promoter ratio")
p1
