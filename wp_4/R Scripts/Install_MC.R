# install bioconductor 
install.packages("BiocManager", dependencies = TRUE)
BiocManager::install(version = "3.14")

# install a few bioconductor dependencies that aren't automatically installed
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'), 
                     dependencies = TRUE, 
                     update = FALSE, 
                     force = TRUE)

# install monocle3 through the cole-trapnell-lab Github, execute
install.packages("devtools", dependencies =  TRUE)
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

# testing the installation
library(monocle3)

# Cicero also requires following packages from bioconductor
#install.packages("BiocManager")
BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"))

# install Cicero
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")

# testing the installation
library(cicero)

