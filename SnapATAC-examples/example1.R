
library(SnapATAC);
x.sp = createSnap(
  file="/home/rstudio/stud3/data/snapATAC_examples/example_1/atac_v1_adult_brain_fresh_5k.snap",
  sample="atac_v1_adult_brain_fresh_5k",
  num.cores=1
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
