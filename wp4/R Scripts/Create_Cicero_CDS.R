library(ggplot2)
library(dplyr)
library(monocle3)
library(cicero)

cds_object <- detect_genes(cds_object)
cds_object <- estimate_size_factors(cds_object)
cds_object <- preprocess_cds(cds_object, 
                             method = "LSI")

cds_object <- reduce_dimension(cds_object, 
                               reduction_method = 'UMAP', 
                               preprocess_method = "LSI", 
                               umap.min_dist = 0.4,
                               umap.n_neighbors=15)

cds_object <- cluster_cells(cds_object)
plot_cells(cds_object, 
           color_cells_by = 'cluster',
           group_label_size = 4,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_groups_by_cluster = TRUE)

umap_coords <- reducedDims(cds_object)$UMAP
head(umap_coords)
cicero_cds <- make_cicero_cds(cds_object, reduced_coordinates = umap_coords)