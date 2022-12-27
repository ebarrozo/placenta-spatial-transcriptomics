# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

set.seed(seed=1)
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
library(SeuratDisk)

######################## ######################## ########################
######################## Had issues installing monocle3 on aagaardlab3, so saved csf-macrophage-subset and doing monocle3 on local
setwd("/Users/enricobarrozo/Library/CloudStorage/Box-Box/AagaardLab/Visium/SpaceRanger/results/seurat_human_v2/subset.macs/")
	## start with cluster 9 based on M0 markers


library(ggplot2)
library(Seurat)
load("human_spatial_data-macs-annotated_v1.RData")


dir.create("Monocle3-lineage")
setwd("Monocle3-lineage")
    ## Also run with "Polarization"
    ## dir.create("Monocle3-Polarization")
    ## setwd("Monocle3-Polarization")

Idents(seurat.object) <- "Lineage"
    ## Idents(seurat.object) <- "Polarization"

dir.create("partition.false")
setwd("partition.false")


library(ggplot2)
library(Seurat)

library(devtools)

library(monocle)
library(monocle3)
library(dplyr)
library(magrittr)
library(data.table)
library(xlsx)

set.seed(seed=1)

length(Seurat::Cells(seurat.object))
        ##13357 ;; 896 amniocyte only subset

##colnames(data)
pData <- seurat.object@meta.data
pData$cell <- rownames(pData)
##pd <- new('AnnotatedDataFrame', data = pData)

##Extract data, phenotype data, and feature data from the SeuratObject
monocle.mat <- as(as.matrix(seurat.object@assays$Spatial@counts), 'sparseMatrix')

##data <- as.matrix(seuratX@assays$integrated@data)
rownames(monocle.mat)
colnames(monocle.mat)
length(colnames(monocle.mat))
        ##13357

##fData <- monocle3::fData(data)
fData <- data.frame(id = row.names(monocle.mat), 
                    gene_short_name = row.names(monocle.mat))
##fd <- new('AnnotatedDataFrame', data = fData)
rownames(fData) <- rownames(monocle.mat)

##Construct monocle cds
seurat.object.monocle <- monocle3::new_cell_data_set(expression_data=monocle.mat,
                                                  cell_metadata = pData,
                                                  gene_metadata = fData)
ncol(seurat.object.monocle)
unique(seurat.object.monocle@colData$Lineage)
seurat.object.cds <- seurat.object.monocle
seurat.object.cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- seurat.object@reductions[["umap"]]@cell.embeddings
seurat.object.cds <- monocle3::cluster_cells(seurat.object.cds,
                                                                                reduction_method="UMAP", 
                                                                                cluster_method="louvain",
                                        verbose=TRUE, 
                                        use_partition=FALSE)
plot_cells(seurat.object.cds)

seurat.object.cds <- monocle3::learn_graph(seurat.object.cds, 
                                        verbose=TRUE, 
                                        use_partition=FALSE)
###Error: No cell clusters for UMAP calculated. Please run cluster_cells with reduction_method = UMAP before running learn_graph.

seurat.object.cds <- monocle3::order_cells(seurat.object.cds, 
                                        reduction_method = "UMAP",
                                        verbose=TRUE)
############################################################################################################################################################
####################################### MANUALLY SELECT START POINT HERE #######################################
############################################################################################################################################################
monocle3_plot_seurat.object_seurat_orig_ident <- monocle3::plot_cells(seurat.object.cds,
                                                                   color_cells_by = "Lineage",
                                                                   label_cell_groups=TRUE,
                                                                   label_leaves=TRUE,
                                                                   label_branch_points=TRUE,
                                                                   graph_label_size=3) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_orig_ident

ggsave("monocle3_plot_seurat.object_seurat_orig_ident_partition_false.pdf", plot = monocle3_plot_seurat.object_seurat_orig_ident, width = 7, height = 5, units = "in")

monocle3_plot_seurat.object_seurat_orig_ident <- monocle3::plot_cells(seurat.object.cds,
                                                                   color_cells_by = "Polarization",
                                                                   label_cell_groups=TRUE,
                                                                   label_leaves=TRUE,
                                                                   label_branch_points=TRUE,
                                                                   graph_label_size=3) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_orig_ident

ggsave("monocle3_plot_seurat.object_seurat_Polarization_partition_false.pdf", plot = monocle3_plot_seurat.object_seurat_orig_ident, width = 7, height = 5, units = "in")

##Rerun using partition = true
##Deleting all of Mike's ggsave commands and replacing them with mine
#ggsave("monocle3_plot_seurat.object_seurat_orig_ident_partition_true.pdf", plot = monocle3_plot_seurat.object_seurat_orig_ident, width = 7, height = 5, units = "in")



monocle3_plot_seurat.object_seurat_partition <- monocle3::plot_cells(seurat.object.cds,
                                                                  color_cells_by = "partition",
                                                                  label_cell_groups=TRUE,
                                                                  label_leaves=TRUE,
                                                                  label_branch_points=TRUE,
                                                                  graph_label_size=3) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_partition

ggsave("monocle_seurat.object_seurat_trajectory_partition_false.pdf", plot = monocle3_plot_seurat.object_seurat_partition, width = 7, height = 5, units = "in")



monocle3_plot_seurat.object_seurat_seurat_clusters <- monocle3::plot_cells(seurat.object.cds,
                                                                        color_cells_by = "seurat_clusters",
                                                                        ##label_cell_groups=TRUE,
                                                                        label_leaves=TRUE,
                                                                        label_branch_points=TRUE,
                                                                        graph_label_size=3,
                                                                        ##group_label_size=5,
                                                                        cell_size = 1) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_seurat_clusters

ggsave("monocle_seurat.object_seurat_trajectory_seurat_clusters_partition_false.pdf", plot = monocle3_plot_seurat.object_seurat_seurat_clusters, width = 7, height = 5, units = "in")


monocle3_plot_seurat.object_seurat_pseudotime <- monocle3::plot_cells(seurat.object.cds,
                                                                   color_cells_by = "pseudotime",
                                                                   label_cell_groups=FALSE,
                                                                   label_leaves=TRUE,
                                                                   label_branch_points=TRUE,
                                                                   graph_label_size=3,
                                                                   cell_size = 1) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_pseudotime
ggsave("monocle_seurat.object_seurat_trajectory_pseudotime_partition_false.pdf", plot = monocle3_plot_seurat.object_seurat_pseudotime, width = 7, height = 5, units = "in")

setwd("..")

####################################### RERUN USING PARTITION = TRUE #######################################

dir.create("partition.true")
setwd("partition.true")

seurat.object.cds <- monocle3::cluster_cells(seurat.object.cds,
                                                                                reduction_method="UMAP", 
                                                                                cluster_method="louvain",
                                        verbose=TRUE, 
                                        use_partition=TRUE)
#plot_cells(seurat.object.cds)                                        

seurat.object.cds <- monocle3::learn_graph(seurat.object.cds, 
                                        verbose=TRUE, 
                                        use_partition=TRUE)

seurat.object.cds <- monocle3::order_cells(seurat.object.cds, 
                                        reduction_method = "UMAP",
                                        verbose=TRUE)
##root_cells remove and select cluster manually with lowest sum of viral mean counts
############################################################################################################################################################
####################################### MANUALLY SELECT START POINT HERE #######################################
############################################################################################################################################################

monocle3_plot_seurat.object_seurat_orig_ident <- monocle3::plot_cells(seurat.object.cds,
                                                                   color_cells_by = "Lineage",
                                                                   label_cell_groups=TRUE,
                                                                   label_leaves=TRUE,
                                                                   label_branch_points=TRUE,
                                                                   graph_label_size=3) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_orig_ident

ggsave("monocle3_plot_seurat.object_seurat_orig_ident_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_orig_ident, width = 7, height = 5, units = "in")

monocle3_plot_seurat.object_seurat_orig_ident <- monocle3::plot_cells(seurat.object.cds,
                                                                   color_cells_by = "Polarization",
                                                                   label_cell_groups=TRUE,
                                                                   label_leaves=TRUE,
                                                                   label_branch_points=TRUE,
                                                                   graph_label_size=3) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_orig_ident

ggsave("monocle3_plot_seurat.object_seurat_Polarization_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_orig_ident, width = 7, height = 5, units = "in")


monocle3_plot_seurat.object_seurat_partition <- monocle3::plot_cells(seurat.object.cds,
                                                                  color_cells_by = "partition",
                                                                  label_cell_groups=TRUE,
                                                                  label_leaves=TRUE,
                                                                  label_branch_points=TRUE,
                                                                  graph_label_size=3) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_partition

ggsave("monocle_seurat.object_seurat_trajectory_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_partition, width = 7, height = 5, units = "in")



monocle3_plot_seurat.object_seurat_seurat_clusters <- monocle3::plot_cells(seurat.object.cds,
                                                                        color_cells_by = "seurat_clusters",
                                                                        ##label_cell_groups=TRUE,
                                                                        label_leaves=TRUE,
                                                                        label_branch_points=TRUE,
                                                                        graph_label_size=3,
                                                                        ##group_label_size=5,
                                                                        cell_size = 1) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_seurat_clusters

ggsave("monocle_seurat.object_seurat_trajectory_seurat_clusters_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_seurat_clusters, width = 7, height = 5, units = "in")

monocle3_plot_seurat.object_seurat_seurat_clusters <- monocle3::plot_cells(seurat.object.cds,
                                                                        color_cells_by = "Lineage",
                                                                        ##label_cell_groups=TRUE,
                                                                        label_leaves=TRUE,
                                                                        label_branch_points=TRUE,
                                                                        graph_label_size=3,
                                                                        ##group_label_size=5,
                                                                        cell_size = 1) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_seurat_clusters

ggsave("monocle_seurat.object_seurat_trajectory_Lineage_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_seurat_clusters, width = 7, height = 5, units = "in")

monocle3_plot_seurat.object_seurat_seurat_clusters <- monocle3::plot_cells(seurat.object.cds,
                                                                        color_cells_by = "Lineage",
                                                                        ##label_cell_groups=TRUE,
                                                                        label_leaves=TRUE,
                                                                        label_branch_points=TRUE,
                                                                        graph_label_size=3,
                                                                        ##group_label_size=5,
                                                                        cell_size = 1) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_seurat_clusters

ggsave("monocle_seurat.object_seurat_trajectory_Lineage_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_seurat_clusters, width = 7, height = 5, units = "in")



monocle3_plot_seurat.object_seurat_pseudotime <- monocle3::plot_cells(seurat.object.cds,
                                                                   color_cells_by = "pseudotime",
                                                                   label_cell_groups=FALSE,
                                                                   label_leaves=TRUE,
                                                                   label_branch_points=TRUE,
                                                                   graph_label_size=3,
                                                                   cell_size = 1) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_pseudotime
ggsave("monocle_seurat.object_seurat_trajectory_pseudotime_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_pseudotime, width = 7, height = 5, units = "in")

setwd("..")
####################################### RERUN USING MONOCLE DIMENSION REDUCTION #######################################
dir.create("monocle.dimred")
setwd("monocle.dimred")


seurat.object.cds <- monocle3::preprocess_cds(seurat.object.cds, num_dim = 50)
#seurat.object.cds <- monocle3::align_cds(seurat.object.cds, preprocess_method = 'PCA', alignment_group = "orig.ident")

seurat.object.cds <- monocle3::reduce_dimension(seurat.object.cds,
                                                                                        preprocess_method = 'PCA')
plot_cells(seurat.object.cds) 

seurat.object.cds <- monocle3::cluster_cells(seurat.object.cds,
                                                                                reduction_method="UMAP", 
                                                                                cluster_method="louvain",
                                        verbose=TRUE, 
                                        use_partition=TRUE)

plot_cells(seurat.object.cds, color_cells_by = "seurat_clusters")


seurat.object.cds <- monocle3::learn_graph(seurat.object.cds, 
                                        verbose=TRUE, 
                                        use_partition=TRUE)

seurat.object.cds <- monocle3::order_cells(seurat.object.cds, 
                                        reduction_method = "UMAP",
                                        verbose=TRUE)
############################################################################################################################################################
####################################### MANUALLY SELECT START POINT HERE #######################################
############################################################################################################################################################


monocle3_plot_seurat.object_seurat_orig_ident <- monocle3::plot_cells(seurat.object.cds,
                                                                   color_cells_by = "Lineage",
                                                                   label_cell_groups=TRUE,
                                                                   label_leaves=TRUE,
                                                                   label_branch_points=TRUE,
                                                                   graph_label_size=3) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_orig_ident

ggsave("monocle3_plot_seurat.object_seurat_orig_ident_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_orig_ident, width = 7, height = 5, units = "in")

monocle3_plot_seurat.object_seurat_partition <- monocle3::plot_cells(seurat.object.cds,
                                                                  color_cells_by = "partition",
                                                                  label_cell_groups=TRUE,
                                                                  label_leaves=TRUE,
                                                                  label_branch_points=TRUE,
                                                                  graph_label_size=3) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_partition

ggsave("monocle_seurat.object_seurat_trajectory_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_partition, width = 7, height = 5, units = "in")



monocle3_plot_seurat.object_seurat_seurat_clusters <- monocle3::plot_cells(seurat.object.cds,
                                                                        color_cells_by = "seurat_clusters",
                                                                        ##label_cell_groups=TRUE,
                                                                        label_leaves=TRUE,
                                                                        label_branch_points=TRUE,
                                                                        graph_label_size=3,
                                                                        ##group_label_size=5,
                                                                        cell_size = 1) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_seurat_clusters

ggsave("monocle_seurat.object_seurat_trajectory_seurat_clusters_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_seurat_clusters, width = 7, height = 5, units = "in")


monocle3_plot_seurat.object_seurat_pseudotime <- monocle3::plot_cells(seurat.object.cds,
                                                                   color_cells_by = "pseudotime",
                                                                   label_cell_groups=FALSE,
                                                                   label_leaves=TRUE,
                                                                   label_branch_points=TRUE,
                                                                   graph_label_size=3,
                                                                   cell_size = 1) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_pseudotime
ggsave("monocle_seurat.object_seurat_trajectory_pseudotime_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_pseudotime, width = 7, height = 5, units = "in")
setwd("..")

####################################### CLUSTERING IN 3D#######################################
dir.create("tweeedeeee")
setwd("tweeedeeee")



########Changing the names of clusters so cluster 0 can be plotted in 3d
##needs to be updated to match the number of clusters you want/have in Seurat
seurat.object <- seurat.object

relative.rank <- seurat.object
new.cluster.ids <- c("z","1","2","3","4")
names(new.cluster.ids) <- levels(relative.rank)
relative.rank <- RenameIdents(relative.rank, new.cluster.ids)

###seurat.object has to be clustered cells
seurat.object <- relative.rank

##Script from Mike that adds seurat_clusters as a way to color your UMAP
monocle.mat <- as(as.matrix(seurat.object@assays$Spatial@counts), 'sparseMatrix')
rownames(monocle.mat)
colnames(monocle.mat)
length(colnames(monocle.mat))
pData <- seurat.object@meta.data
pData$cell <- rownames(pData)
fData <- data.frame(id = row.names(monocle.mat), 
                    gene_short_name = row.names(monocle.mat))
rownames(fData) <- rownames(monocle.mat)
seurat.object.monocle <- monocle3::new_cell_data_set(expression_data=monocle.mat,
                                                  cell_metadata = pData,
                                                  gene_metadata = fData)
ncol(seurat.object.monocle)
unique(seurat.object.monocle@colData$orig.ident)
seurat.object.cds <- seurat.object.monocle
seurat.object.cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- seurat.object@reductions[["umap"]]@cell.embeddings

#clustering in 3d
seurat.object.cds <- monocle3::preprocess_cds(seurat.object.cds, num_dim = 50)
seurat.object.cds.3d <- seurat.object.cds
seurat.object.cds.3d <- reduce_dimension(seurat.object.cds, max_components = 3)
seurat.object.cds.3d <- cluster_cells(seurat.object.cds.3d)
plot_cells_3d(seurat.object.cds.3d, color_cells_by="Lineage")
####################################### FIND YOUR START POINT BEFORE PROCEEDING
######### export as webpage-clusters
seurat.object.cds.3d <- learn_graph(seurat.object.cds.3d)
############################################################################################################################################################
####################################### MANUALLY SELECT START POINT HERE #######################################
############################################################################################################################################################
###### I'm just picking the tail... 
seurat.object.cds.3d <- order_cells(seurat.object.cds.3d)

d.seurat <- plot_cells_3d(seurat.object.cds.3d, color_cells_by="Lineage")
d.seurat
#Save .html manually
#ggsave("3d.seurat.pdf", plot = d.seurat, width = 7, height = 5, units = "in")

d.pseudotime <- plot_cells_3d(seurat.object.cds.3d, color_cells_by="pseudotime")
d.pseudotime
#Save .html manually
#ggsave("3d_pseudotime.pdf", plot = d.seurat, width = 7, height = 5, units = "in")

####################################### GO HOME AND HAVE A DRINK #######################################
setwd("..")


############## 9.24.20_ MM- More Pruned ###############################################################
##defaults() list(1,1/3,10,FALSE,TRUE)
graph.control.params <- list(1,1/3,20,FALSE,TRUE)
names(graph.control.params) <- c("euclidean_distance_ratio",
                                 "geodesic_distance_ratio",
                                 "minimal_branch_len",
                                 "orthogonal_proj_tip",
                                 "prune_graph")
seurat.object.cds.3d.more.pruned <- learn_graph(seurat.object.cds.3d,
                                       use_partition = TRUE,
                                       close_loop = TRUE,
                                       learn_graph_control = graph.control.params ##default is 10
                            )
seurat.object.cds.3d.more.pruned <- order_cells(seurat.object.cds.3d.more.pruned)
d.seurat.more.pruned <- plot_cells_3d(seurat.object.cds.3d.more.pruned,
                                 color_cells_by="seurat_clusters")
d.seurat.more.pruned
####################################### PRUNING #######################################
####################################### PRUNING #######################################
####################################### PRUNING #######################################
## viral-optimized: highest viral counts 7>1>10>8>5>2>6>0>3>4>13>12>9>15>11>14>16
## Start at cluster 7;; 16
dir.create("partition.true-pruned-20")
setwd("partition.true-pruned-20")
##Extract data, phenotype data, and feature data from the SeuratObject
monocle.mat <- as(as.matrix(seurat.object@assays$Spatial@counts), 'sparseMatrix')
##data <- as.matrix(seuratX@assays$integrated@data)
rownames(monocle.mat)
colnames(monocle.mat)
length(colnames(monocle.mat))
##5535

##colnames(data)
pData <- seurat.object@meta.data
pData$cell <- rownames(pData)
##pd <- new('AnnotatedDataFrame', data = pData)

##fData <- monocle3::fData(data)
fData <- data.frame(id = row.names(monocle.mat), 
                    gene_short_name = row.names(monocle.mat))
##fd <- new('AnnotatedDataFrame', data = fData)
rownames(fData) <- rownames(monocle.mat)

##Construct monocle cds
seurat.object.monocle <- monocle3::new_cell_data_set(expression_data=monocle.mat,
                                                  cell_metadata = pData,
                                                  gene_metadata = fData)
##lowerDetectionLimit = 0.5,
##expressionFamily = uninormal()) ## since I have already normalized, thresholded and scaled in Suerat v3.0.0.9150

ncol(seurat.object.monocle)
unique(seurat.object.monocle@colData$orig.ident)

seurat.object.cds <- seurat.object.monocle

seurat.object.cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- seurat.object@reductions[["umap"]]@cell.embeddings
#seurat.object.cds@clusters$listData[["UMAP"]] <- seurat.object@reductions[["umap"]]@cell.embeddings

seurat.object.cds <- monocle3::cluster_cells(seurat.object.cds,
                                                                                reduction_method="UMAP", 
                                                                                cluster_method="louvain",
                                        verbose=TRUE, 
                                        use_partition=TRUE)
#plot_cells(seurat.object.cds)                                        
####################################### THIS IS WHERE PRUNING HAPPENS

graph.control.params <- list(1,1/3,20,FALSE,TRUE)
names(graph.control.params) <- c("euclidean_distance_ratio",
                                 "geodesic_distance_ratio",
                                 "minimal_branch_len",
                                 "orthogonal_proj_tip",
                                 "prune_graph")

seurat.object.cds <- monocle3::learn_graph(seurat.object.cds, 
                                        use_partition = TRUE,
                                       close_loop = TRUE,
                                       learn_graph_control = graph.control.params ##default is 10
                            )

seurat.object.cds <- monocle3::order_cells(seurat.object.cds, 
                                        reduction_method = "UMAP",
                                        verbose=TRUE)
##root_cells remove and select cluster manually with lowest sum of viral mean counts
############################################################################################################################################################
####################################### MANUALLY SELECT START POINT HERE #######################################
############################################################################################################################################################

####################################### IT TAKES A MINUTE, BE PATIENT

monocle3_plot_seurat.object_seurat_orig_ident <- monocle3::plot_cells(seurat.object.cds,
                                                                   color_cells_by = "Lineage",
                                                                   label_cell_groups=TRUE,
                                                                   label_leaves=TRUE,
                                                                   label_branch_points=TRUE,
                                                                   graph_label_size=3) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_orig_ident

ggsave("monocle3_plot_seurat.object_seurat_Lineage_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_orig_ident, width = 7, height = 5, units = "in")

monocle3_plot_seurat.object_seurat_orig_ident <- monocle3::plot_cells(seurat.object.cds,
                                                                   color_cells_by = "Polarization",
                                                                   label_cell_groups=TRUE,
                                                                   label_leaves=TRUE,
                                                                   label_branch_points=TRUE,
                                                                   graph_label_size=3) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_orig_ident

ggsave("monocle3_plot_seurat.object_seurat_Polarization_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_orig_ident, width = 7, height = 5, units = "in")


monocle3_plot_seurat.object_seurat_partition <- monocle3::plot_cells(seurat.object.cds,
                                                                  color_cells_by = "partition",
                                                                  label_cell_groups=TRUE,
                                                                  label_leaves=TRUE,
                                                                  label_branch_points=TRUE,
                                                                  graph_label_size=3) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_partition

ggsave("monocle_seurat.object_seurat_trajectory_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_partition, width = 7, height = 5, units = "in")



monocle3_plot_seurat.object_seurat_seurat_clusters <- monocle3::plot_cells(seurat.object.cds,
                                                                        color_cells_by = "seurat_clusters",
                                                                        ##label_cell_groups=TRUE,
                                                                        label_leaves=TRUE,
                                                                        label_branch_points=TRUE,
                                                                        graph_label_size=3,
                                                                        ##group_label_size=5,
                                                                        cell_size = 1) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_seurat_clusters

ggsave("monocle_seurat.object_seurat_trajectory_seurat_clusters_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_seurat_clusters, width = 7, height = 5, units = "in")


monocle3_plot_seurat.object_seurat_pseudotime <- monocle3::plot_cells(seurat.object.cds,
                                                                   color_cells_by = "pseudotime",
                                                                   label_cell_groups=FALSE,
                                                                   label_leaves=TRUE,
                                                                   label_branch_points=TRUE,
                                                                   graph_label_size=3,
                                                                   cell_size = 1) +
  theme(legend.position = "right")

monocle3_plot_seurat.object_seurat_pseudotime
ggsave("monocle_seurat.object_seurat_trajectory_pseudotime_partition_TRUE.pdf", plot = monocle3_plot_seurat.object_seurat_pseudotime, width = 7, height = 5, units = "in")

setwd("..")

####################################################################################################