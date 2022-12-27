## visium_human_seurat_v2.R

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

    ## v1 was samples 01, 02, 03 and 15. V2 includes those samples and S16, S17, S18, S19, S20, S21, S22, S23a, S23b, S24, S25, and S26. 
    # Metadata
      ## uninfected: S01, S02, S03, S24, S25, S26
      ## SARS-CoV-2 asymptomatic: S15, S17, S18, S19
      ## SARS-CoV-2 symptomatic: S16, S20, S21, S22, S23a, S23b
        ## SARS-CoV-2 perivillous fibrinoids: S23a and S23b
        ## Noted areas of inflammation: S18
      ## Samples S01, S02, and S03 were from the sample placenta from distinct regions including chorionic villi, decidua, and a chorioamiotic membrane roll
      ## Samples 15-26 were sampled from the placental parenchyme with various regions of villi and some decidua
      ## Samples 23a and 23b were from the same placenta punch, selected for two sections due to the presence of perivillious fibrinoid pathology identified in QC H&E slide

## Analyze data on AagaardLab3
# 

set.seed(seed=1)
setwd("/home/ebarrozo/visium/results")
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
library(SeuratDisk)
library(dplyr)
library(patchwork)
dir.create("seurat_human_v2")
setwd("seurat_human_v2")

setwd("/home/ebarrozo/visium/results/seurat_human_v2")

############################################################# ############################################################# #############################################################
############################################################# S01: Load data, Quality Control (QC) filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S01-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S01.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S01", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S01.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S01.manual <- RenameGenesSeurat(obj = S01.manual, newnames = new.names)
row.names(S01.manual)

## Save a table of all.genes
all.genes <- row.names(S01.manual)
write.table(all.genes, "S01.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S01.manual <- PercentageFeatureSet(S01.manual, pattern = "^MT-", col.name = "percent.mt")
S01.manual <- PercentageFeatureSet(S01.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S01.manual <- PercentageFeatureSet(S01.manual, features = SARS.genes, col.name = "percent.viral")
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S01.manual<- NormalizeData(S01.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S01.manual <- CellCycleScoring(S01.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S01.manual$CC.Difference <- S01.manual$S.Score - S01.manual$G2M.Score

acute.viral <- subset(S01.manual, features = rownames(S01.manual)[rownames(S01.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S01.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 548
acute.viral <- subset(acute.viral, percent.viral > 0)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(S01.manual)[rownames(S01.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

S01.manual$orig.ident <- "S01"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S01.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S01.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S01.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S01.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S01.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S01.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S01.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S01.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S01.manual)
  #  548
S01.manual2 <- subset(S01.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S01.manual2)
  #   511
qc.vlnplot <- VlnPlot(S01.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S01.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S01.manual <- S01.manual2
rm(S01.manual2)
qc.vlnplot <- VlnPlot(S01.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S01.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S01.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S01.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S01.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S01: Dimension Reduction, Visualization, and Differential Expression (DE) #############################################################
############################################################# ############################################################# #############################################################

## Run SCTransform to normalize by neg. binomial
S01.manual <- SCTransform(S01.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S01.manual <- RunPCA(S01.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S01.manual)
ggsave("S01.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S01.manual <- FindNeighbors(S01.manual, reduction = "pca", dims = 1:30)
S01.manual <- FindClusters(S01.manual, verbose = T)
S01.manual <- RunUMAP(S01.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S01.manual, reduction = "umap", label = TRUE)
p1
  # clusters 0-5
p2 <- SpatialDimPlot(S01.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S01.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
  # As there are many colors, it can be challenging to visualize which voxel belongs to which cluster. We have a few strategies to help with this. Setting the label parameter places a colored box at the median of each cluster (see the plot above).
p3 <- wrap_plots(p1, p2, p4)
ggsave("S01.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S01.manual, cells.highlight = CellsByIdentities(object = S01.manual, idents = c(0, 1, 2, 3, 4, 5)), facet.highlight = TRUE, ncol = 3)
ggsave("S01.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S01.manual <- FindSpatiallyVariableFeatures(S01.manual, assay = "SCT", features = VariableFeatures(S01.manual)[1:1000], selection.method = "markvariogram")
S01.manual.var <- S01.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S01.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S01.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S01.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S01.manual <- NormalizeData(S01.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S01.manual <- ScaleData(S01.manual, features = S01.manual.var, assay="Spatial")

Idents(S01.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S01.manual, features = S01.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S01_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
top5.heatmap <- DoHeatmap(S01.manual, features = top5$gene, raster = FALSE)
ggsave("S01.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S01.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S01.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S01.manual, features = unique.top2)
png(file=paste0("S01.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S01.manual, features = unique.top2)
png(file=paste0("S01.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S01.manual, features = unique.top2)
png(file=paste0("S01.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)


############################################################# ############################################################# #############################################################
############################################################# S03: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################
## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S03-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S03.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S03", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S03.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S03.manual <- RenameGenesSeurat(obj = S03.manual, newnames = new.names)
row.names(S03.manual)

## Save a table of all.genes
all.genes <- row.names(S03.manual)
write.table(all.genes, "S03.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S03.manual <- PercentageFeatureSet(S03.manual, pattern = "^MT-", col.name = "percent.mt")
S03.manual <- PercentageFeatureSet(S03.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S03.manual <- PercentageFeatureSet(S03.manual, features = SARS.genes, col.name = "percent.viral")
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S03.manual<- NormalizeData(S03.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S03.manual <- CellCycleScoring(S03.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S03.manual$CC.Difference <- S03.manual$S.Score - S03.manual$G2M.Score

acute.viral <- subset(S03.manual, features = rownames(S03.manual)[rownames(S03.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S03.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 470
acute.viral <- subset(acute.viral, percent.viral > 0)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(S03.manual)[rownames(S03.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

S03.manual$orig.ident <- "S03"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S03.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S03.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S03.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S03.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S03.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S03.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S03.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S03.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S03.manual)
  #  470
S03.manual2 <- subset(S03.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S03.manual2)
  #   455
qc.vlnplot <- VlnPlot(S03.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S03.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S03.manual <- S03.manual2
rm(S03.manual2)
qc.vlnplot <- VlnPlot(S03.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S03.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S03.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S03.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S03.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S03: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S03.manual <- SCTransform(S03.manual, assay = "Spatial", verbose = FALSE)

## Dimensionality reduction, clustering, and visualization
S03.manual <- RunPCA(S03.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S03.manual)
ggsave("S03.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S03.manual <- FindNeighbors(S03.manual, reduction = "pca", dims = 1:30)
S03.manual <- FindClusters(S03.manual, verbose = T)
S03.manual <- RunUMAP(S03.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S03.manual, reduction = "umap", label = TRUE)
p1
  # clusters 0-6
p2 <- SpatialDimPlot(S03.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S03.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S03.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
  ## Change idents to match the number of clusters
p3<- SpatialDimPlot(S03.manual, cells.highlight = CellsByIdentities(object = S03.manual, idents = c(0, 1, 2, 3, 4, 5, 6)), facet.highlight = TRUE, ncol = 3)
ggsave("S03.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S03.manual <- FindSpatiallyVariableFeatures(S03.manual, assay = "SCT", features = VariableFeatures(S03.manual)[1:1000], selection.method = "markvariogram")
S03.manual.var <- S03.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S03.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S03.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S03.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S03.manual <- NormalizeData(S03.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S03.manual <- ScaleData(S03.manual, features = S03.manual.var, assay="Spatial")

Idents(S03.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S03.manual, features = S03.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S03_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
top5.heatmap <- DoHeatmap(S03.manual, features = top5$gene, raster = FALSE)
ggsave("S03.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S03.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S03.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S03.manual, features = unique.top2)
png(file=paste0("S03.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S03.manual, features = unique.top2)
png(file=paste0("S03.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S03.manual, features = unique.top2)
png(file=paste0("S03.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# S04: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################
## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S04-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S04.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S04", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S04.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S04.manual <- RenameGenesSeurat(obj = S04.manual, newnames = new.names)
row.names(S04.manual)

## Save a table of all.genes
all.genes <- row.names(S04.manual)
write.table(all.genes, "S04.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S04.manual <- PercentageFeatureSet(S04.manual, pattern = "^MT-", col.name = "percent.mt")
S04.manual <- PercentageFeatureSet(S04.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S04.manual <- PercentageFeatureSet(S04.manual, features = SARS.genes, col.name = "percent.viral")
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S04.manual<- NormalizeData(S04.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S04.manual <- CellCycleScoring(S04.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S04.manual$CC.Difference <- S04.manual$S.Score - S04.manual$G2M.Score

acute.viral <- subset(S04.manual, features = rownames(S04.manual)[rownames(S04.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S04.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 548
acute.viral <- subset(acute.viral, percent.viral > 0)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(S04.manual)[rownames(S04.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

S04.manual$orig.ident <- "S04"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S04.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S04.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S04.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S04.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S04.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S04.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S04.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S04.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S04.manual)
  #  967
S04.manual2 <- subset(S04.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S04.manual2)
  #   951
qc.vlnplot <- VlnPlot(S04.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S04.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S04.manual <- S04.manual2
rm(S04.manual2)
qc.vlnplot <- VlnPlot(S04.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S04.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S04.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S04.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S04.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S04: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S04.manual <- SCTransform(S04.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S04.manual <- RunPCA(S04.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S04.manual)
ggsave("S04.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S04.manual <- FindNeighbors(S04.manual, reduction = "pca", dims = 1:30)
S04.manual <- FindClusters(S04.manual, verbose = T)
S04.manual <- RunUMAP(S04.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S04.manual, reduction = "umap", label = TRUE)
p1
  # clusters 0-6
p2 <- SpatialDimPlot(S04.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S04.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S04.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S04.manual, cells.highlight = CellsByIdentities(object = S04.manual, idents = c(0, 1, 2, 3, 4, 5, 6)), facet.highlight = TRUE, ncol = 3)
ggsave("S04.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S04.manual <- FindSpatiallyVariableFeatures(S04.manual, assay = "SCT", features = VariableFeatures(S04.manual)[1:1000], selection.method = "markvariogram")
S04.manual.var <- S04.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S04.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S04.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S04.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S04.manual <- NormalizeData(S04.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S04.manual <- ScaleData(S04.manual, features = S04.manual.var, assay="Spatial")

Idents(S04.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S04.manual, features = S04.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S04_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
top5.heatmap <- DoHeatmap(S04.manual, features = top5$gene, raster = FALSE)
ggsave("S04.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S04.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S04.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S04.manual, features = unique.top2)
png(file=paste0("S04.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S04.manual, features = unique.top2)
png(file=paste0("S04.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S04.manual, features = unique.top2)
png(file=paste0("S04.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# S15: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S15-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S15.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S15", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S15.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S15.manual <- RenameGenesSeurat(obj = S15.manual, newnames = new.names)
row.names(S15.manual)

## Save a table of all.genes
all.genes <- row.names(S15.manual)
write.table(all.genes, "S15.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S15.manual <- PercentageFeatureSet(S15.manual, pattern = "^MT-", col.name = "percent.mt")
S15.manual <- PercentageFeatureSet(S15.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S15.manual <- PercentageFeatureSet(S15.manual, features = SARS.genes, col.name = "percent.viral")
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S15.manual<- NormalizeData(S15.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S15.manual <- CellCycleScoring(S15.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S15.manual$CC.Difference <- S15.manual$S.Score - S15.manual$G2M.Score

acute.viral <- subset(S15.manual, features = rownames(S15.manual)[rownames(S15.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S15.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 1505
acute.viral <- subset(acute.viral, percent.viral > 0)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(S15.manual)[rownames(S15.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

S15.manual$orig.ident <- "S15"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S15.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S15.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S15.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S15.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S15.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S15.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S15.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S15.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S15.manual)
  #  1505
S15.manual2 <- subset(S15.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S15.manual2)
  #   1480
qc.vlnplot <- VlnPlot(S15.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S15.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S15.manual <- S15.manual2
rm(S15.manual2)
qc.vlnplot <- VlnPlot(S15.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S15.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S15.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S15.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S15.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S15: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S15.manual <- SCTransform(S15.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S15.manual <- RunPCA(S15.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S15.manual)
ggsave("S15.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S15.manual <- FindNeighbors(S15.manual, reduction = "pca", dims = 1:30)
S15.manual <- FindClusters(S15.manual, verbose = T)
S15.manual <- RunUMAP(S15.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S15.manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(S15.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S15.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S15.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S15.manual, cells.highlight = CellsByIdentities(object = S15.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave("S15.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S15.manual <- FindSpatiallyVariableFeatures(S15.manual, assay = "SCT", features = VariableFeatures(S15.manual)[1:1000], selection.method = "markvariogram")
S15.manual.var <- S15.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S15.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S15.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S15.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S15.manual <- NormalizeData(S15.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S15.manual <- ScaleData(S15.manual, features = S15.manual.var, assay="Spatial")

Idents(S15.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S15.manual, features = S15.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S15_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S15.manual, features = S15.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S15.manual, features = top5$gene, raster = FALSE)
ggsave("S15.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S15.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S15.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S15.manual, features = unique.top2)
png(file=paste0("S15.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S15.manual, features = unique.top2)
png(file=paste0("S15.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S15.manual, features = unique.top2)
png(file=paste0("S15.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# S16: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S16-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S16.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S16", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S16.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S16.manual <- RenameGenesSeurat(obj = S16.manual, newnames = new.names)
row.names(S16.manual)

## Save a table of all.genes
all.genes <- row.names(S16.manual)
write.table(all.genes, "S16.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S16.manual <- PercentageFeatureSet(S16.manual, pattern = "^MT-", col.name = "percent.mt")
S16.manual <- PercentageFeatureSet(S16.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S16.manual <- PercentageFeatureSet(S16.manual, features = SARS.genes, col.name = "percent.viral")
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S16.manual<- NormalizeData(S16.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S16.manual <- CellCycleScoring(S16.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S16.manual$CC.Difference <- S16.manual$S.Score - S16.manual$G2M.Score

acute.viral <- subset(S16.manual, features = rownames(S16.manual)[rownames(S16.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S16.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 1505
acute.viral <- subset(acute.viral, percent.viral > 0)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(S16.manual)[rownames(S16.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

S16.manual$orig.ident <- "S16"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S16.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S16.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S16.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S16.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S16.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S16.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S16.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S16.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S16.manual)
  #  1547
S16.manual2 <- subset(S16.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S16.manual2)
  #   1512
qc.vlnplot <- VlnPlot(S16.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S16.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S16.manual <- S16.manual2
rm(S16.manual2)
qc.vlnplot <- VlnPlot(S16.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S16.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S16.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S16.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S16.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S16: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S16.manual <- SCTransform(S16.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S16.manual <- RunPCA(S16.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S16.manual)
ggsave("S16.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S16.manual <- FindNeighbors(S16.manual, reduction = "pca", dims = 1:30)
S16.manual <- FindClusters(S16.manual, verbose = T)
S16.manual <- RunUMAP(S16.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S16.manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(S16.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S16.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S16.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S16.manual, cells.highlight = CellsByIdentities(object = S16.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave("S16.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S16.manual <- FindSpatiallyVariableFeatures(S16.manual, assay = "SCT", features = VariableFeatures(S16.manual)[1:1000], selection.method = "markvariogram")
S16.manual.var <- S16.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S16.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S16.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S16.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S16.manual <- NormalizeData(S16.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S16.manual <- ScaleData(S16.manual, features = S16.manual.var, assay="Spatial")

Idents(S16.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S16.manual, features = S16.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S16_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S16.manual, features = S16.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S16.manual, features = top5$gene, raster = FALSE)
ggsave("S16.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S16.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S16.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S16.manual, features = unique.top2)
png(file=paste0("S16.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S16.manual, features = unique.top2)
png(file=paste0("S16.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S16.manual, features = unique.top2)
png(file=paste0("S16.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# S17: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S17-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S17.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S17", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S17.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S17.manual <- RenameGenesSeurat(obj = S17.manual, newnames = new.names)
row.names(S17.manual)

## Save a table of all.genes
all.genes <- row.names(S17.manual)
write.table(all.genes, "S17.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S17.manual <- PercentageFeatureSet(S17.manual, pattern = "^MT-", col.name = "percent.mt")
S17.manual <- PercentageFeatureSet(S17.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S17.manual <- PercentageFeatureSet(S17.manual, features = SARS.genes, col.name = "percent.viral")
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S17.manual<- NormalizeData(S17.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S17.manual <- CellCycleScoring(S17.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S17.manual$CC.Difference <- S17.manual$S.Score - S17.manual$G2M.Score

acute.viral <- subset(S17.manual, features = rownames(S17.manual)[rownames(S17.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S17.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 1375
acute.viral <- subset(acute.viral, percent.viral > 0)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(S17.manual)[rownames(S17.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

S17.manual$orig.ident <- "S17"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S17.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S17.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S17.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S17.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S17.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S17.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S17.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S17.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S17.manual)
  #  1375
S17.manual2 <- subset(S17.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S17.manual2)
  #   1356
qc.vlnplot <- VlnPlot(S17.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S17.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S17.manual <- S17.manual2
rm(S17.manual2)
qc.vlnplot <- VlnPlot(S17.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S17.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S17.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S17.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S17.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S17: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S17.manual <- SCTransform(S17.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S17.manual <- RunPCA(S17.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S17.manual)
ggsave("S17.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S17.manual <- FindNeighbors(S17.manual, reduction = "pca", dims = 1:30)
S17.manual <- FindClusters(S17.manual, verbose = T)
S17.manual <- RunUMAP(S17.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S17.manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(S17.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S17.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S17.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S17.manual, cells.highlight = CellsByIdentities(object = S17.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave("S17.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S17.manual <- FindSpatiallyVariableFeatures(S17.manual, assay = "SCT", features = VariableFeatures(S17.manual)[1:1000], selection.method = "markvariogram")
S17.manual.var <- S17.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S17.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S17.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S17.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S17.manual <- NormalizeData(S17.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S17.manual <- ScaleData(S17.manual, features = S17.manual.var, assay="Spatial")

Idents(S17.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S17.manual, features = S17.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S17_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S17.manual, features = S17.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S17.manual, features = top5$gene, raster = FALSE)
ggsave("S17.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S17.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S17.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S17.manual, features = unique.top2)
png(file=paste0("S17.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S17.manual, features = unique.top2)
png(file=paste0("S17.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S17.manual, features = unique.top2)
png(file=paste0("S17.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# S18: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S18-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S18.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S18", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S18.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S18.manual <- RenameGenesSeurat(obj = S18.manual, newnames = new.names)
row.names(S18.manual)

## Save a table of all.genes
all.genes <- row.names(S18.manual)
write.table(all.genes, "S18.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S18.manual <- PercentageFeatureSet(S18.manual, pattern = "^MT-", col.name = "percent.mt")
S18.manual <- PercentageFeatureSet(S18.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S18.manual <- PercentageFeatureSet(S18.manual, features = SARS.genes, col.name = "percent.viral")
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S18.manual<- NormalizeData(S18.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S18.manual <- CellCycleScoring(S18.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S18.manual$CC.Difference <- S18.manual$S.Score - S18.manual$G2M.Score

acute.viral <- subset(S18.manual, features = rownames(S18.manual)[rownames(S18.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S18.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 1561
acute.viral <- subset(acute.viral, percent.viral > 0)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(S18.manual)[rownames(S18.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

S18.manual$orig.ident <- "S18"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S18.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S18.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S18.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S18.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S18.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S18.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S18.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S18.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
qc.vlnplot

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S18.manual)
  #  1561
S18.manual2 <- subset(S18.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S18.manual2)
  #   1498
qc.vlnplot <- VlnPlot(S18.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S18.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S18.manual <- S18.manual2
rm(S18.manual2)
qc.vlnplot <- VlnPlot(S18.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S18.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S18.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S18.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S18.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S18: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S18.manual <- SCTransform(S18.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S18.manual <- RunPCA(S18.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S18.manual)
ggsave("S18.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S18.manual <- FindNeighbors(S18.manual, reduction = "pca", dims = 1:30)
S18.manual <- FindClusters(S18.manual, verbose = T)
S18.manual <- RunUMAP(S18.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S18.manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(S18.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S18.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S18.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S18.manual, cells.highlight = CellsByIdentities(object = S18.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave("S18.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S18.manual <- FindSpatiallyVariableFeatures(S18.manual, assay = "SCT", features = VariableFeatures(S18.manual)[1:1000], selection.method = "markvariogram")
S18.manual.var <- S18.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S18.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S18.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S18.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S18.manual <- NormalizeData(S18.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S18.manual <- ScaleData(S18.manual, features = S18.manual.var, assay="Spatial")

Idents(S18.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S18.manual, features = S18.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S18_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S18.manual, features = S18.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S18.manual, features = top5$gene, raster = FALSE)
ggsave("S18.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S18.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S18.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S18.manual, features = unique.top2)
png(file=paste0("S18.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S18.manual, features = unique.top2)
png(file=paste0("S18.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S18.manual, features = unique.top2)
png(file=paste0("S18.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# S19: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S19-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S19.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S19", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S19.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S19.manual <- RenameGenesSeurat(obj = S19.manual, newnames = new.names)
row.names(S19.manual)

## Save a table of all.genes
all.genes <- row.names(S19.manual)
write.table(all.genes, "S19.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S19.manual <- PercentageFeatureSet(S19.manual, pattern = "^MT-", col.name = "percent.mt")
S19.manual <- PercentageFeatureSet(S19.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S19.manual <- PercentageFeatureSet(S19.manual, features = SARS.genes, col.name = "percent.viral")
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S19.manual<- NormalizeData(S19.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S19.manual <- CellCycleScoring(S19.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S19.manual$CC.Difference <- S19.manual$S.Score - S19.manual$G2M.Score

acute.viral <- subset(S19.manual, features = rownames(S19.manual)[rownames(S19.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S19.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 918
acute.viral <- subset(acute.viral, percent.viral > 0)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(S19.manual)[rownames(S19.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

S19.manual$orig.ident <- "S19"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S19.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S19.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S19.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S19.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S19.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S19.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S19.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S19.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S19.manual)
  #  918
S19.manual2 <- subset(S19.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S19.manual2)
  #   910
qc.vlnplot <- VlnPlot(S19.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S19.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S19.manual <- S19.manual2
rm(S19.manual2)
qc.vlnplot <- VlnPlot(S19.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S19.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S19.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S19.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S19.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S19: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S19.manual <- SCTransform(S19.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S19.manual <- RunPCA(S19.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S19.manual)
ggsave("S19.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S19.manual <- FindNeighbors(S19.manual, reduction = "pca", dims = 1:30)
S19.manual <- FindClusters(S19.manual, verbose = T)
S19.manual <- RunUMAP(S19.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S19.manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(S19.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S19.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S19.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S19.manual, cells.highlight = CellsByIdentities(object = S19.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave("S19.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S19.manual <- FindSpatiallyVariableFeatures(S19.manual, assay = "SCT", features = VariableFeatures(S19.manual)[1:1000], selection.method = "markvariogram")
S19.manual.var <- S19.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S19.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S19.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S19.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S19.manual <- NormalizeData(S19.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S19.manual <- ScaleData(S19.manual, features = S19.manual.var, assay="Spatial")

Idents(S19.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S19.manual, features = S19.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S19_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S19.manual, features = S19.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S19.manual, features = top5$gene, raster = FALSE)
ggsave("S19.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S19.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S19.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S19.manual, features = unique.top2)
png(file=paste0("S19.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S19.manual, features = unique.top2)
png(file=paste0("S19.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S19.manual, features = unique.top2)
png(file=paste0("S19.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# S20: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S20-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S20.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S20", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S20.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S20.manual <- RenameGenesSeurat(obj = S20.manual, newnames = new.names)
row.names(S20.manual)

## Save a table of all.genes
all.genes <- row.names(S20.manual)
write.table(all.genes, "S20.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S20.manual <- PercentageFeatureSet(S20.manual, pattern = "^MT-", col.name = "percent.mt")
S20.manual <- PercentageFeatureSet(S20.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S20.manual <- PercentageFeatureSet(S20.manual, features = SARS.genes, col.name = "percent.viral")
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S20.manual<- NormalizeData(S20.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S20.manual <- CellCycleScoring(S20.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S20.manual$CC.Difference <- S20.manual$S.Score - S20.manual$G2M.Score

acute.viral <- subset(S20.manual, features = rownames(S20.manual)[rownames(S20.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S20.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 515
acute.viral <- subset(acute.viral, percent.viral > 0)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(S20.manual)[rownames(S20.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

S20.manual$orig.ident <- "S20"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S20.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S20.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S20.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S20.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S20.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S20.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S20.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S20.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S20.manual)
  #  515
S20.manual2 <- subset(S20.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S20.manual2)
  #   492
qc.vlnplot <- VlnPlot(S20.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S20.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S20.manual <- S20.manual2
rm(S20.manual2)
qc.vlnplot <- VlnPlot(S20.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S20.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S20.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S20.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S20.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S20: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S20.manual <- SCTransform(S20.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S20.manual <- RunPCA(S20.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S20.manual)
ggsave("S20.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S20.manual <- FindNeighbors(S20.manual, reduction = "pca", dims = 1:30)
S20.manual <- FindClusters(S20.manual, verbose = T)
S20.manual <- RunUMAP(S20.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S20.manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(S20.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S20.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S20.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S20.manual, cells.highlight = CellsByIdentities(object = S20.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave("S20.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S20.manual <- FindSpatiallyVariableFeatures(S20.manual, assay = "SCT", features = VariableFeatures(S20.manual)[1:1000], selection.method = "markvariogram")
S20.manual.var <- S20.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S20.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S20.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S20.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S20.manual <- NormalizeData(S20.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S20.manual <- ScaleData(S20.manual, features = S20.manual.var, assay="Spatial")

Idents(S20.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S20.manual, features = S20.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S20_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S20.manual, features = S20.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S20.manual, features = top5$gene, raster = FALSE)
ggsave("S20.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S20.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S20.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S20.manual, features = unique.top2)
png(file=paste0("S20.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S20.manual, features = unique.top2)
png(file=paste0("S20.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S20.manual, features = unique.top2)
png(file=paste0("S20.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# S21: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S21-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S21.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S21", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S21.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S21.manual <- RenameGenesSeurat(obj = S21.manual, newnames = new.names)
row.names(S21.manual)

## Save a table of all.genes
all.genes <- row.names(S21.manual)
write.table(all.genes, "S21.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S21.manual <- PercentageFeatureSet(S21.manual, pattern = "^MT-", col.name = "percent.mt")
S21.manual <- PercentageFeatureSet(S21.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S21.manual <- PercentageFeatureSet(S21.manual, features = SARS.genes, col.name = "percent.viral")
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S21.manual<- NormalizeData(S21.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S21.manual <- CellCycleScoring(S21.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S21.manual$CC.Difference <- S21.manual$S.Score - S21.manual$G2M.Score

acute.viral <- subset(S21.manual, features = rownames(S21.manual)[rownames(S21.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S21.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 341
acute.viral <- subset(acute.viral, percent.viral > 0)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(S21.manual)[rownames(S21.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

S21.manual$orig.ident <- "S21"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S21.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S21.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S21.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S21.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S21.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S21.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S21.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S21.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S21.manual)
  #  341
S21.manual2 <- subset(S21.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S21.manual2)
  #   326
qc.vlnplot <- VlnPlot(S21.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S21.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S21.manual <- S21.manual2
rm(S21.manual2)
qc.vlnplot <- VlnPlot(S21.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S21.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S21.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S21.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S21.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S21: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S21.manual <- SCTransform(S21.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S21.manual <- RunPCA(S21.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S21.manual)
ggsave("S21.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S21.manual <- FindNeighbors(S21.manual, reduction = "pca", dims = 1:30)
S21.manual <- FindClusters(S21.manual, verbose = T)
S21.manual <- RunUMAP(S21.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S21.manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(S21.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S21.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S21.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S21.manual, cells.highlight = CellsByIdentities(object = S21.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave("S21.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S21.manual <- FindSpatiallyVariableFeatures(S21.manual, assay = "SCT", features = VariableFeatures(S21.manual)[1:1000], selection.method = "markvariogram")
S21.manual.var <- S21.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S21.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S21.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S21.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S21.manual <- NormalizeData(S21.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S21.manual <- ScaleData(S21.manual, features = S21.manual.var, assay="Spatial")

Idents(S21.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S21.manual, features = S21.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S21_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S21.manual, features = S21.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S21.manual, features = top5$gene, raster = FALSE)
ggsave("S21.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S21.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S21.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S21.manual, features = unique.top2)
png(file=paste0("S21.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S21.manual, features = unique.top2)
png(file=paste0("S21.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S21.manual, features = unique.top2)
png(file=paste0("S21.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# S22: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S22-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S22.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S22", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S22.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S22.manual <- RenameGenesSeurat(obj = S22.manual, newnames = new.names)
row.names(S22.manual)

## Save a table of all.genes
all.genes <- row.names(S22.manual)
write.table(all.genes, "S22.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S22.manual <- PercentageFeatureSet(S22.manual, pattern = "^MT-", col.name = "percent.mt")
S22.manual <- PercentageFeatureSet(S22.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S22.manual <- PercentageFeatureSet(S22.manual, features = SARS.genes, col.name = "percent.viral")
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S22.manual<- NormalizeData(S22.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S22.manual <- CellCycleScoring(S22.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S22.manual$CC.Difference <- S22.manual$S.Score - S22.manual$G2M.Score

acute.viral <- subset(S22.manual, features = rownames(S22.manual)[rownames(S22.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S22.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 669
acute.viral <- subset(acute.viral, percent.viral > 0)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(S22.manual)[rownames(S22.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

S22.manual$orig.ident <- "S22"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S22.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S22.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S22.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S22.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S22.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S22.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S22.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S22.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S22.manual)
  #  669
S22.manual2 <- subset(S22.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S22.manual2)
  #   661
qc.vlnplot <- VlnPlot(S22.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S22.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S22.manual <- S22.manual2
rm(S22.manual2)
qc.vlnplot <- VlnPlot(S22.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S22.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S22.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S22.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S22.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S22: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S22.manual <- SCTransform(S22.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S22.manual <- RunPCA(S22.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S22.manual)
ggsave("S22.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S22.manual <- FindNeighbors(S22.manual, reduction = "pca", dims = 1:30)
S22.manual <- FindClusters(S22.manual, verbose = T)
S22.manual <- RunUMAP(S22.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S22.manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(S22.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S22.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S22.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S22.manual, cells.highlight = CellsByIdentities(object = S22.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave("S22.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S22.manual <- FindSpatiallyVariableFeatures(S22.manual, assay = "SCT", features = VariableFeatures(S22.manual)[1:1000], selection.method = "markvariogram")
S22.manual.var <- S22.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S22.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S22.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S22.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S22.manual <- NormalizeData(S22.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S22.manual <- ScaleData(S22.manual, features = S22.manual.var, assay="Spatial")

Idents(S22.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S22.manual, features = S22.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S22_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S22.manual, features = S22.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S22.manual, features = top5$gene, raster = FALSE)
ggsave("S22.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S22.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S22.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S22.manual, features = unique.top2)
png(file=paste0("S22.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S22.manual, features = unique.top2)
png(file=paste0("S22.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S22.manual, features = unique.top2)
png(file=paste0("S22.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

######################################################### ############################################################# #############################################################
############################################################# S23a: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S23a-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S23a.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S23a", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S23a.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S23a.manual <- RenameGenesSeurat(obj = S23a.manual, newnames = new.names)
row.names(S23a.manual)

## Save a table of all.genes
all.genes <- row.names(S23a.manual)
write.table(all.genes, "S23a.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S23a.manual <- PercentageFeatureSet(S23a.manual, pattern = "^MT-", col.name = "percent.mt")
S23a.manual <- PercentageFeatureSet(S23a.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S23a.manual <- PercentageFeatureSet(S23a.manual, features = SARS.genes, col.name = "percent.viral")


## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S23a.manual<- NormalizeData(S23a.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S23a.manual <- CellCycleScoring(S23a.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S23a.manual$CC.Difference <- S23a.manual$S.Score - S23a.manual$G2M.Score

acute.viral <- subset(S23a.manual, features = rownames(S23a.manual)[rownames(S23a.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S23a.viral.counts_all-cells.txt", sep="\t")

ncol(acute.viral)   # 369
acute.viral <- subset(acute.viral, percent.viral > 0)
ncol(acute.viral)   # 15
acute.viral <- subset(acute.viral, features = rownames(S23a.manual)[rownames(S23a.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S23a.viral.counts_infected-cells.txt", sep="\t")

S23a.manual$orig.ident <- "S23a"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S23a.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S23a.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S23a.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S23a.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S23a.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S23a.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S23a.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S23a.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S23a.manual)
  #  369

S23a.manual2 <- subset(S23a.manual, subset = nFeature_Spatial > 5 & nFeature_Spatial < 10000 & nCount_Spatial > 50 & nCount_Spatial < 100000  & percent.mt < 30 & percent.ribo < 20)
ncol(S23a.manual2)
  #   93 - changing filters to match this sample # S23a.manual2 <- subset(S23a.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
  #   260 - filters changed above to reflect quality changes due to suspected viral infection. subset = nFeature_Spatial > 100 & nFeature_Spatial < 10000 & nCount_Spatial > 200 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 20)
    ##  361 - using filters from 23b
qc.vlnplot <- VlnPlot(S23a.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S23a.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S23a.manual <- S23a.manual2
rm(S23a.manual2)
qc.vlnplot <- VlnPlot(S23a.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S23a.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S23a.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S23a.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S23a.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S23a: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S23a.manual <- SCTransform(S23a.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S23a.manual <- RunPCA(S23a.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S23a.manual)
ggsave("S23a.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S23a.manual <- FindNeighbors(S23a.manual, reduction = "pca", dims = 1:30)
S23a.manual <- FindClusters(S23a.manual, verbose = T)
S23a.manual <- RunUMAP(S23a.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S23a.manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(S23a.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S23a.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S23a.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S23a.manual, cells.highlight = CellsByIdentities(object = S23a.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave("S23a.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S23a.manual <- FindSpatiallyVariableFeatures(S23a.manual, assay = "SCT", features = VariableFeatures(S23a.manual)[1:1000], selection.method = "markvariogram")
S23a.manual.var <- S23a.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S23a.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S23a.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S23a.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S23a.manual <- NormalizeData(S23a.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S23a.manual <- ScaleData(S23a.manual, features = S23a.manual.var, assay="Spatial")

Idents(S23a.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S23a.manual, features = S23a.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S23a_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S23a.manual, features = S23a.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S23a.manual, features = top5$gene, raster = FALSE)
ggsave("S23a.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S23a.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S23a.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S23a.manual, features = unique.top2)
png(file=paste0("S23a.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S23a.manual, features = unique.top2)
png(file=paste0("S23a.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S23a.manual, features = unique.top2)
png(file=paste0("S23a.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# S23b: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S23b-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S23b.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S23b", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S23b.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S23b.manual <- RenameGenesSeurat(obj = S23b.manual, newnames = new.names)
row.names(S23b.manual)

## Save a table of all.genes
all.genes <- row.names(S23b.manual)
write.table(all.genes, "S23b.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S23b.manual <- PercentageFeatureSet(S23b.manual, pattern = "^MT-", col.name = "percent.mt")
S23b.manual <- PercentageFeatureSet(S23b.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S23b.manual <- PercentageFeatureSet(S23b.manual, features = SARS.genes, col.name = "percent.viral", assay = 'Spatial')
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing
S23b.manual <- PercentageFeatureSet(S23b.manual, features =c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10"), col.name = "percent.viral", assay = 'Spatial')


## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S23b.manual<- NormalizeData(S23b.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S23b.manual <- CellCycleScoring(S23b.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S23b.manual$CC.Difference <- S23b.manual$S.Score - S23b.manual$G2M.Score

acute.viral <- subset(S23b.manual, features = rownames(S23b.manual)[rownames(S23b.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S23b.viral.counts_all-cells.txt", sep="\t")
ncol(acute.viral)   # 1323 total cells
acute.viral <- subset(acute.viral, percent.viral > 0)
ncol(acute.viral)   # 1068 cells with SARS-CoV-2 transcripts
acute.viral <- subset(acute.viral, features = rownames(S23b.manual)[rownames(S23b.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S23b.viral.counts_infected-cells.txt", sep="\t")

S23b.manual$orig.ident <- "S23b"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S23b.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S23b.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S23b.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S23b.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S23b.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S23b.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S23b.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S23b.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S23b.manual)
  #  1323
S23b.manual2 <- subset(S23b.manual, subset = nFeature_Spatial > 5 & nFeature_Spatial < 10000 & nCount_Spatial > 50 & nCount_Spatial < 100000  & percent.mt < 30 & percent.ribo < 20)
ncol(S23b.manual2)
  #    - only 1 cell left with current filters. subset(S23b.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
    ## 850-  Because these are very infected cells, adjusted filters. Most shockingly, the nFeatures average is <100 and average UMIs/cell is well below 500

qc.vlnplot <- VlnPlot(S23b.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S23b.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S23b.manual <- S23b.manual2
rm(S23b.manual2)
qc.vlnplot <- VlnPlot(S23b.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S23b.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S23b.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S23b.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S23b.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S23b: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S23b.manual <- SCTransform(S23b.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S23b.manual <- RunPCA(S23b.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S23b.manual)
ggsave("S23b.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S23b.manual <- FindNeighbors(S23b.manual, reduction = "pca", dims = 1:30)
S23b.manual <- FindClusters(S23b.manual, verbose = T)
S23b.manual <- RunUMAP(S23b.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S23b.manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(S23b.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S23b.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S23b.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S23b.manual, cells.highlight = CellsByIdentities(object = S23b.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave("S23b.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S23b.manual <- FindSpatiallyVariableFeatures(S23b.manual, assay = "SCT", features = VariableFeatures(S23b.manual)[1:1000], selection.method = "markvariogram")
S23b.manual.var <- S23b.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S23b.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S23b.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S23b.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S23b.manual <- NormalizeData(S23b.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S23b.manual <- ScaleData(S23b.manual, features = S23b.manual.var, assay="Spatial")

Idents(S23b.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S23b.manual, features = S23b.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S23b_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S23b.manual, features = S23b.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S23b.manual, features = top5$gene, raster = FALSE)
ggsave("S23b.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S23b.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S23b.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S23b.manual, features = unique.top2)
png(file=paste0("S23b.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S23b.manual, features = unique.top2)
png(file=paste0("S23b.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S23b.manual, features = unique.top2)
png(file=paste0("S23b.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# S24: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S24-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S24.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S24", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S24.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S24.manual <- RenameGenesSeurat(obj = S24.manual, newnames = new.names)
row.names(S24.manual)

## Save a table of all.genes
all.genes <- row.names(S24.manual)
write.table(all.genes, "S24.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S24.manual <- PercentageFeatureSet(S24.manual, pattern = "^MT-", col.name = "percent.mt")
S24.manual <- PercentageFeatureSet(S24.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S24.manual <- PercentageFeatureSet(S24.manual, features = SARS.genes, col.name = "percent.viral")


## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S24.manual<- NormalizeData(S24.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S24.manual <- CellCycleScoring(S24.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S24.manual$CC.Difference <- S24.manual$S.Score - S24.manual$G2M.Score

acute.viral <- subset(S24.manual, features = rownames(S24.manual)[rownames(S24.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S24.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 1505
acute.viral <- subset(acute.viral, percent.viral > 0.0000000001)
  # Error: No cells found
acute.viral <- subset(acute.viral, features = rownames(S24.manual)[rownames(S24.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S24.viral.counts_infected-cells.txt", sep="\t")

S24.manual$orig.ident <- "S24"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S24.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S24.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S24.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S24.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S24.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S24.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S24.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S24.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S24.manual)
  #  2027
S24.manual2 <- subset(S24.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S24.manual2)
  #   1756
qc.vlnplot <- VlnPlot(S24.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S24.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S24.manual <- S24.manual2
rm(S24.manual2)
qc.vlnplot <- VlnPlot(S24.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S24.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S24.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S24.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S24.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S24: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S24.manual <- SCTransform(S24.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S24.manual <- RunPCA(S24.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S24.manual)
ggsave("S24.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S24.manual <- FindNeighbors(S24.manual, reduction = "pca", dims = 1:30)
S24.manual <- FindClusters(S24.manual, verbose = T)
S24.manual <- RunUMAP(S24.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S24.manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(S24.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S24.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S24.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S24.manual, cells.highlight = CellsByIdentities(object = S24.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave("S24.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S24.manual <- FindSpatiallyVariableFeatures(S24.manual, assay = "SCT", features = VariableFeatures(S24.manual)[1:1000], selection.method = "markvariogram")
S24.manual.var <- S24.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S24.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S24.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S24.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S24.manual <- NormalizeData(S24.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S24.manual <- ScaleData(S24.manual, features = S24.manual.var, assay="Spatial")

Idents(S24.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S24.manual, features = S24.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S24_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S24.manual, features = S24.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S24.manual, features = top5$gene, raster = FALSE)
ggsave("S24.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S24.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S24.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S24.manual, features = unique.top2)
png(file=paste0("S24.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S24.manual, features = unique.top2)
png(file=paste0("S24.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S24.manual, features = unique.top2)
png(file=paste0("S24.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# S25: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S25-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S25.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S25", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S25.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S25.manual <- RenameGenesSeurat(obj = S25.manual, newnames = new.names)
row.names(S25.manual)

## Save a table of all.genes
all.genes <- row.names(S25.manual)
write.table(all.genes, "S25.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S25.manual <- PercentageFeatureSet(S25.manual, pattern = "^MT-", col.name = "percent.mt")
S25.manual <- PercentageFeatureSet(S25.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S25.manual <- PercentageFeatureSet(S25.manual, features = SARS.genes, col.name = "percent.viral")

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S25.manual<- NormalizeData(S25.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S25.manual <- CellCycleScoring(S25.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S25.manual$CC.Difference <- S25.manual$S.Score - S25.manual$G2M.Score

acute.viral <- subset(S25.manual, features = rownames(S25.manual)[rownames(S25.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S25.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 1505
acute.viral <- subset(acute.viral, percent.viral > 0)
# Error: No cells found
ncol(acute.viral)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(S25.manual)[rownames(S25.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

S25.manual$orig.ident <- "S25"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S25.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S25.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S25.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S25.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S25.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S25.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S25.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S25.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S25.manual)
  #  2605
S25.manual2 <- subset(S25.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S25.manual2)
  #   2488
qc.vlnplot <- VlnPlot(S25.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S25.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S25.manual <- S25.manual2
rm(S25.manual2)
qc.vlnplot <- VlnPlot(S25.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S25.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S25.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S25.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S25.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S25: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S25.manual <- SCTransform(S25.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S25.manual <- RunPCA(S25.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S25.manual)
ggsave("S25.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S25.manual <- FindNeighbors(S25.manual, reduction = "pca", dims = 1:30)
S25.manual <- FindClusters(S25.manual, verbose = T)
S25.manual <- RunUMAP(S25.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S25.manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(S25.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S25.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S25.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S25.manual, cells.highlight = CellsByIdentities(object = S25.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave("S25.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S25.manual <- FindSpatiallyVariableFeatures(S25.manual, assay = "SCT", features = VariableFeatures(S25.manual)[1:1000], selection.method = "markvariogram")
S25.manual.var <- S25.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S25.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S25.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S25.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S25.manual <- NormalizeData(S25.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S25.manual <- ScaleData(S25.manual, features = S25.manual.var, assay="Spatial")

Idents(S25.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S25.manual, features = S25.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S25_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S25.manual, features = S25.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S25.manual, features = top5$gene, raster = FALSE)
ggsave("S25.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S25.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S25.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S25.manual, features = unique.top2)
png(file=paste0("S25.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S25.manual, features = unique.top2)
png(file=paste0("S25.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S25.manual, features = unique.top2)
png(file=paste0("S25.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

############################################################# ############################################################# #############################################################
############################################################# S26: Load data, QC filtering #############################################################
############################################################# ############################################################# #############################################################

## https://satijalab.org/seurat/reference/index.html#section-spatial
data_dir <- '/home/ebarrozo/visium/results/S26-manual/outs'
list.files(data_dir)
#### scripts modified from https://satijalab.org/seurat/articles/spatial_vignette.html#10x-visium  compiled: 2021-08-30
S26.manual <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "S26", filter.matrix = TRUE, to.upper = TRUE)

## Fix transcript names that have extra dashes added after custom ref generation (e.g. GRCh38----------------------)
old.names <- row.names(S26.manual)
new.names <- gsub ("GRCH38------(.*)", "\\1", old.names)
new.names <- gsub ("NC-045512.2-(.*)", "\\1", new.names)
new.names <- gsub ("GENE-(.*)", "\\1", new.names)

# Function from https://github.com/satijalab/seurat/issues/2617; modified for Spatial assay instead of RNA
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = new.names) { 
  print("Run this before integration. It only changes obj@assays$Spatial@counts, @data and @scale.data.")
  Spatial <- obj@assays$Spatial
  if (nrow(Spatial) == length(newnames)) {
    if (length(Spatial@counts)) Spatial@counts@Dimnames[[1]]            <- newnames
    if (length(Spatial@data)) Spatial@data@Dimnames[[1]]                <- newnames
    if (length(Spatial@scale.data)) Spatial@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(Spatial) != nrow(newnames)"}
  obj@assays$Spatial <- Spatial
  return(obj)
}

S26.manual <- RenameGenesSeurat(obj = S26.manual, newnames = new.names)
row.names(S26.manual)

## Save a table of all.genes
all.genes <- row.names(S26.manual)
write.table(all.genes, "S26.all.genes.txt", sep="\t")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

## Assign and examine unfiltered quality control (QC) metrics
S26.manual <- PercentageFeatureSet(S26.manual, pattern = "^MT-", col.name = "percent.mt")
S26.manual <- PercentageFeatureSet(S26.manual, pattern = "^RP[SL]", col.name = "percent.ribo")
S26.manual <- PercentageFeatureSet(S26.manual, features = SARS.genes, col.name = "percent.viral")
  # Error in h(simpleError(msg, call)) : 
    # error in evaluating the argument 'x' in selecting a method for function 'colSums': invalid character indexing

## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
S26.manual<- NormalizeData(S26.manual, normalization.method = "LogNormalize", scale.factor = 10000)


# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene

# Run cell cycle scoring 
S26.manual <- CellCycleScoring(S26.manual, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S26.manual$CC.Difference <- S26.manual$S.Score - S26.manual$G2M.Score

acute.viral <- subset(S26.manual, features = rownames(S26.manual)[rownames(S26.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "S26.viral.counts_all-cells.txt", sep="\t")

  #### Confirmed manually there are no SARS-CoV-2 transcripts ###########  =SUM(B32287:PN32287)

ncol(acute.viral)   # 1505
acute.viral <- subset(acute.viral, percent.viral > 0)
  # Error in FetchData(object = object, vars = unique(x = expr.char[vars.use]),  : 
    # None of the requested variables were found: 
acute.viral <- subset(acute.viral, features = rownames(S26.manual)[rownames(S26.manual) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")

S26.manual$orig.ident <- "S26"   

## plot UMI counts/spot 
plot1 <- VlnPlot(S26.manual, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S26.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
    # pt.size.factor = 1.6 default
    # alpha - minimum and maximum transparency. Default is c(1, 1).
    # Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
    ## crop = FALSE,
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S26.manual_unfiltered_nCount_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

plot1 <- VlnPlot(S26.manual, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(S26.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S26.manual_unfiltered_nFeature_Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

qc.vlnplot <- VlnPlot(S26.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S26.manual_unfiltered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

qc.vlnplot
## Iteratively change filter settings based on what this library looks like. Examine unfiltered vs filtered to see the final effects (filtering out likely doublets/empty GEMs)
ncol(S26.manual)
  #  1505
S26.manual2 <- subset(S26.manual, subset = nFeature_Spatial > 200 & nFeature_Spatial < 10000 & nCount_Spatial > 1000 & nCount_Spatial < 100000 & percent.mt > 0.000001 & percent.mt < 25 & percent.ribo < 15)
ncol(S26.manual2)
  #   1480
qc.vlnplot <- VlnPlot(S26.manual2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S26.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")
S26.manual <- S26.manual2
rm(S26.manual2)
qc.vlnplot <- VlnPlot(S26.manual, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), ncol=4, group.by = "orig.ident", pt.size = 0.000001) + NoLegend()
ggsave("S26.manual_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 16, height = 8, units = "in")

plot1 <- SpatialFeaturePlot(S26.manual, features = "nCount_Spatial", pt.size.factor = 2.0) + theme(legend.position = "right")
plot2 <- SpatialFeaturePlot(S26.manual, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
plot3 <- wrap_plots(plot1, plot2)
ggsave("S26.manual_filtered__Spatial.pdf", plot = plot3, device = "pdf", width = 8, height = 4, units = "in")

############################################################# ############################################################# #############################################################
############################################################# S26: Dimension Reduction, Visualization, and DE #############################################################
############################################################# ############################################################# #############################################################
## Run SCTransform to normalize by neg. binomial
S26.manual <- SCTransform(S26.manual, assay = "Spatial", verbose = FALSE)


## Dimensionality reduction, clustering, and visualization
S26.manual <- RunPCA(S26.manual, assay = "SCT", verbose = T)

elbow.plot <- ElbowPlot(S26.manual)
ggsave("S26.manual_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

S26.manual <- FindNeighbors(S26.manual, reduction = "pca", dims = 1:30)
S26.manual <- FindClusters(S26.manual, verbose = T)
S26.manual <- RunUMAP(S26.manual, reduction = "pca", dims = 1:30)
p1 <- DimPlot(S26.manual, reduction = "umap", label = TRUE)
p1
  ##### clusters 0 thru 7
p2 <- SpatialDimPlot(S26.manual, label = TRUE, label.size = 3)
p2
p4 <- SpatialDimPlot(S26.manual, label = TRUE, label.size = 5, pt.size.factor = 4.0)
p4
p3 <- wrap_plots(p1, p2, p4)
ggsave("S26.manual_UMAP.pdf", plot = p3, device = "pdf", width = 12, height = 4, units = "in")

## UMAP split.by clusters
p3<- SpatialDimPlot(S26.manual, cells.highlight = CellsByIdentities(object = S26.manual, idents = c(0, 1, 2, 3, 4, 5, 6, 7)), facet.highlight = TRUE, ncol = 3)
ggsave("S26.manual_UMAP_Spatial-split.by.cluster.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


## Find top variable transcripts. This may take a long time.
S26.manual <- FindSpatiallyVariableFeatures(S26.manual, assay = "SCT", features = VariableFeatures(S26.manual)[1:1000], selection.method = "markvariogram")
S26.manual.var <- S26.manual@assays[["SCT"]]@var.features

top.features <- head(SpatiallyVariableFeatures(S26.manual, selection.method = "markvariogram"), 6)
p3 <- SpatialFeaturePlot(S26.manual, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("S26.manual_spatial.var.features.pdf", plot = p3, device = "pdf", width = 8, height = 4, units = "in")


##### Perform differential expression between clusters
S26.manual <- NormalizeData(S26.manual, normalization.method = "LogNormalize", scale.factor = 10000)
S26.manual <- ScaleData(S26.manual, features = S26.manual.var, assay="Spatial")

Idents(S26.manual) <- "seurat_clusters"
de_markers <- FindAllMarkers(S26.manual, features = S26.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "S26_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
  ## make sure there are some transcripts that passed the lnFC(>0.693)
    ## cluster 0 does not have sig. genes
de_markers <- FindAllMarkers(S26.manual, features = S26.manual.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)

top5.heatmap <- DoHeatmap(S26.manual, features = top5$gene, raster = FALSE)
ggsave("S26.manual_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(S26.manual, features = top5$gene, slot="counts", raster = FALSE)
ggsave("S26.manual_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S26.manual, features = unique.top2)
png(file=paste0("S26.manual_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S26.manual, features = unique.top2)
png(file=paste0("S26.manual_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(S26.manual, features = unique.top2)
png(file=paste0("S26.manual_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)

rm(acute.viral)
rm(plot2)
rm(plot1)
rm(plot3)
rm(qc.vlnplot)
rm(old.names)
rm(new.names)
rm(data_dir)

############################################################# ############################################################# #############################################################
############################################################# Integration: merge all data using CCA anchors #############################################################
############################################################# ############################################################# #############################################################

## Slim down the seurat.objects to save space. Don't need to keep the scale.data
S01.manual2 <- DietSeurat(S01.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
View(S01.manual2)
  ## from 188MB to 166 MB
S01.manual <- S01.manual2
rm(S01.manual2)
S03.manual <- DietSeurat(S03.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S04.manual <- DietSeurat(S04.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S15.manual <- DietSeurat(S15.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S16.manual <- DietSeurat(S16.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S17.manual <- DietSeurat(S17.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S18.manual <- DietSeurat(S18.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S19.manual <- DietSeurat(S19.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S20.manual <- DietSeurat(S20.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S21.manual <- DietSeurat(S21.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S22.manual <- DietSeurat(S22.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S23a.manual <- DietSeurat(S23a.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S23b.manual <- DietSeurat(S23b.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S24.manual <- DietSeurat(S24.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S25.manual <- DietSeurat(S25.manual, counts=TRUE, data=TRUE, scale.data = FALSE)
S26.manual <- DietSeurat(S26.manual, counts=TRUE, data=TRUE, scale.data = FALSE)




## combine lists of top variable features for later DE analysis/clustering
var.combined <- union(S01.manual.var, S03.manual.var)
var.combined <- union(var.combined, S04.manual.var)
var.combined <- union(var.combined, S15.manual.var)
var.combined <- union(var.combined, S16.manual.var)
var.combined <- union(var.combined, S17.manual.var)
var.combined <- union(var.combined, S18.manual.var)
var.combined <- union(var.combined, S19.manual.var)
var.combined <- union(var.combined, S20.manual.var)
var.combined <- union(var.combined, S21.manual.var)
var.combined <- union(var.combined, S22.manual.var)
var.combined <- union(var.combined, S23a.manual.var)
var.combined <- union(var.combined, S23b.manual.var)
var.combined <- union(var.combined, S24.manual.var)
var.combined <- union(var.combined, S25.manual.var)
var.combined <- union(var.combined, S26.manual.var)

  ## 16200 combined variable features

## Merge objects
all_merged <- merge(x = S01.manual, y = c(S01.manual, S03.manual, S04.manual, S15.manual, S16.manual, S17.manual, S18.manual, S19.manual, S20.manual, S21.manual, S22.manual, S23a.manual, S23b.manual, S24.manual, S25.manual, S26.manual), merge.data = TRUE, project = "all_merged")
## make combined QC plot 
qc.vlnplot <- VlnPlot(all_merged, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), group.by = "orig.ident", pt.size = 0.000001, ncol = 4) + NoLegend()
ggsave("all_merged_filtered_QC_vlnplot.pdf", plot = qc.vlnplot, device = "pdf", width = 15, height = 5, units = "in")
rm(qc.vlnplot)
rm(S01.manual)
rm(S03.manual)
rm(S04.manual)
rm(S15.manual)
rm(S01.manual.var)
rm(S03.manual.var)
rm(S04.manual.var)
rm(S15.manual.var)
rm(S16.manual.var)
rm(S17.manual.var)
rm(S18.manual.var)
rm(S19.manual.var)
rm(S20.manual.var)
rm(S21.manual.var)
rm(S22.manual.var)
rm(S23.manual.var)
rm(S24.manual.var)
rm(S25.manual.var)
rm(S26.manual.var)

options(future.globals.maxSize = 300000000000)

all_merged <- SCTransform(all_merged, assay = "Spatial", verbose = T)
all_merged.list <- SplitObject(all_merged, split.by = "orig.ident")


# Make a reference list
reference.list <- all_merged.list[c("S01", "S03", "S04", "S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", "S23a", "S23b", "S24", "S25", "S26")]

# Select transcripts used for data integration
reference.features <- SelectIntegrationFeatures(object.list = reference.list, nfeatures = 3000)
  # No variable features found for object1 in the object.list. Running FindVariableFeatures ...
  # Error: SCT assay is comprised of multiple SCT models. To change the variable features, please set manually with VariableFeatures<-
    ## re-run SCTransform on the all_merged object
write.table(reference.features, "all_merged.reference.features.txt", sep="\t")


reference.list <- PrepSCTIntegration(object.list = reference.list, anchor.features = reference.features, verbose = T)

# Find integration anchors, this takes about 15-20 mins
all_merged.anchors <- FindIntegrationAnchors(object.list = reference.list, normalization.method = "SCT", anchor.features = reference.features, verbose = T)
## If there is an object with few cells/spots (<50?), you will see the error message below and need to change the k.filter to 1 fewer than tne lowest number of cells/spots you have
  # You're computing too large a percentage of total singular values, use a standard svd instead.
  # Warning messages:
    # 1: In FilterAnchors(object = object.pair, assay = assay, slot = slot,  :
      # Number of anchor cells is less than k.filter. Retaining all anchors.
      # all_merged.anchors <- FindIntegrationAnchors(object.list = reference.list, reduction="rpca", normalization.method = "SCT", anchor.features = reference.features, verbose = TRUE)
          ## Projecting new data onto SVD
rm(reference.features)
rm(reference.list)

seurat.object <- IntegrateData(anchorset = all_merged.anchors, normalization.method = "SCT")
  # Merging dataset 3 into 1
  #Extracting anchors for merged samples
  #Finding integration vectors
  #Finding integration vector weights
  #Error in idx[i, ] <- res[[i]][[1]] : 
  #  number of items to replace is not a multiple of replacement length
    #all_merged.integrated <- IntegrateData(anchorset = all_merged.anchors, normalization.method = "SCT", k.weight = 9)
      # Usually the k.weight is set to 100 by default, but here the smallest dataset only has 46 cells. Based on searching with the error message I found this shouldn't be a problem https://github.com/satijalab/seurat/issues/3930 https://github.com/satijalab/seurat/issues/1472
rm("all_merged")
rm(all_merged.list)
rm("all_merged.anchors")

# Add all viral+human genes as a list of features to call upon during clustering or differential expression of all genes, which may be compuationally demanding
all.genes <- rownames(seurat.object@assays[["Spatial"]])

# Use the 'integrated' assay only for clustering
DefaultAssay(seurat.object) <- "integrated"
# Add all integrated genes as a list of features to call upon during clustering (or quick- differential expression)
integrated.genes <- rownames(seurat.object)
write.table(integrated.genes, "integrated.genes.txt", sep="\t")

## Use the 'Spatial' object for QC and differential expression analysis of allllll genes, which is computationally demanding
#DefaultAssay(seurat.object) <- "Spatial"

rm(S16.manual)
rm(S17.manual)
rm(S18.manual)
rm(S19.manual)
rm(S20.manual)
rm(S21.manual)
rm(S22.manual)
rm(S23a.manual)
rm(S23b.manual)
rm(S24.manual)
rm(S25.manual)
rm(S26.manual)
rm(S23a.manual.var)
rm(S23b.manual.var)

setwd("/home/ebarrozo/visium/results/seurat_human_v2")
save.image("human_spatial_data-integrated_v1.RData")
# load("human_spatial_data-integrated_v1.RData")

# Metadata
      ## uninfected: S01, S02, S03, S24, S25, S26
      ## SARS-CoV-2 asymptomatic: S15, S17, S18, S19
      ## SARS-CoV-2 symptomatic: S16, S20, S21, S22, S23a, S23b
        ## SARS-CoV-2 perivillous fibrinoids: S23a and S23b
        ## Noted areas of inflammation: S18
      ## Samples S01, S02, and S03 were from the sample placenta from distinct regions including chorionic villi, decidua, and a chorioamiotic membrane roll
      ## Samples 15-26 were sampled from the placental parenchyme with various regions of villi and some decidua
      ## Samples 23a and 23b were from the same placenta punch, selected for two sections due to the presence of perivillious fibrinoid pathology identified in QC H&E slide


new.metadata <- c("Uninfected", "Uninfected", "Uninfected", "SARS-CoV-2.ax", "SARS-CoV-2.sx", "SARS-CoV-2.ax", "SARS-CoV-2.ax", "SARS-CoV-2.ax", "SARS-CoV-2.sx", "SARS-CoV-2.sx", "SARS-CoV-2.sx", "SARS-CoV-2.sx", "SARS-CoV-2.sx", "Uninfected","Uninfected","Uninfected")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
## add cluster Subtypes as a metadata column
seurat.object$Condition <- Idents(seurat.object)
Idents(seurat.object) <- "Condition"
Idents(seurat.object)


############################################################# ############################################################# #############################################################
############################################################# Integrated: Dimension Reduction and Visualization #############################################################
############################################################# ############################################################# #############################################################

seurat.object <- ScaleData(seurat.object, features = integrated.genes, assay = "integrated", vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"))
seurat.object<- RunPCA(seurat.object, features = integrated.genes, assay = "integrated", slot = "scale.data")
elbow.plot <- ElbowPlot(seurat.object)
ggsave("seurat.object_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

DefaultAssay(seurat.object) <- "integrated"
seurat.object <- FindNeighbors(seurat.object, features = "integrated.genes", dims = 1:30)
seurat.object <- FindClusters(seurat.object, resolution = 1.2)
seurat.object <- RunUMAP(seurat.object, dims=1:30)
  ## clusters 0 thru 10, res=0.6
      ## DE analysis of clusters 0-10 suggested cluster 9 had immune cells that could not be separated from other clusters. Going to try to increase resolution to capture known cellular complexity (an approach used by Vento-Tormo et al.)
        ## res=1.2, clusters 0 thru 15
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE")
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Condition")
p4 <- FeaturePlot(seurat.object, reduction = "umap", features = "percent.viral")

umap.combined <- CombinePlots(plots = list(p1, p2, p3, p4))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("integrated_UMAP.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "orig.ident", label=TRUE)
ggsave("integrated_UMAP_splitby_libID.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "Condition", label=TRUE)
ggsave("integrated_UMAP_splitby_Condition.pdf", plot = p2, device = "pdf", width = 15, height = 3, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "seurat_clusters", label=TRUE)
ggsave("integrated_UMAP_splitby_clusters.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
#examine UMAPs with qc metrics 
p5 <- FeaturePlot(seurat.object, features = 'nFeature_Spatial')
p6 <- FeaturePlot(seurat.object, features = 'nCount_Spatial')
p7 <- FeaturePlot(seurat.object, features = 'percent.mt')
p8 <- FeaturePlot(seurat.object, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")


###### UMAP + Spatial DimPlots
p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)

p2 <- SpatialDimPlot(seurat.object, group.by = "percent.viral", facet.highlight = TRUE, ncol = 3)
p2 <- SpatialDimPlot(seurat.object, features = "percent.viral", ncol = 3)
plot2 <- SpatialFeaturePlot(seurat.object, features = "percent.viral", ncol = 1) + theme(legend.position = "right") + theme(plot.title = element_text(size = rel(1.5), angle = 90))


plot2 <- SpatialFeaturePlot(seurat.object, features = "percent.viral", ncol = 1) + patchwork::plot_layout(ncol = 1)
####### FOR SOME REASON, INTEGRATION ADDS IMAGES ARBITRARILY
ggsave("integrated_SpatialDimPlot_percent.viral_TEST.pdf", plot = plot2, device = "pdf", width = 45, height = 11, units = "in")

ggsave("integrated_SpatialDimPlot_percent.viral_TEST.pdf", plot = plot2, device = "pdf", width = 8, height = 45, units = "in")
seurat.object <- PercentageFeatureSet(seurat.object, features = SARS.genes, col.name = "percent.viral", assay="Spatial")

## Using the test images above, copy and paste the images that work below. 
    ## identified a pattern and confirmed these all worked. 
S01
S03.1
S04.2
S15.3
S16.4
S17.5
S18.6
S19.7
S20.8
S21.9
S22.10
S23a.11
S23b.12
S24.13
S25.14
S26.15

image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_percent.viral.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
plot2<- SpatialDimPlot(seurat.object, images=image.list) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_cropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialDimPlot(seurat.object, images=image.list, crop=FALSE) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")



## Spatial DimPlots split.by clusters
image.list <- c("S01")
p2 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p2
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S01.pdf", plot = p2, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S03.1")
p3 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p3
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S03.pdf", plot = p3, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S04.2")
p4 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p4
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S04.pdf", plot = p4, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S15.3")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p5
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S15.pdf", plot = p5, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S16.4")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p5
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S16.pdf", plot = p5, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S17.5")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p5
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S17.pdf", plot = p5, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S18.6")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p5
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S18.pdf", plot = p5, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S19.7")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p5
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S19.pdf", plot = p5, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S20.8")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p5
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S20.pdf", plot = p5, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S21.9")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p5
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S21.pdf", plot = p5, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S22.10")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p5
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S22.pdf", plot = p5, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S23a.11")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p5
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S23a.pdf", plot = p5, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S23b.12")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p5
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S23b.pdf", plot = p5, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S24.13")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p5
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S24.pdf", plot = p5, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S25.14")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p5
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S25.pdf", plot = p5, device = "pdf", width = 8, height = 15, units = "in")
image.list <- c("S26.15")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)), facet.highlight = TRUE, ncol = 5)
p5
ggsave("ntegrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster-S26.pdf", plot = p5, device = "pdf", width = 8, height = 15, units = "in")


rm(umap.combined)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(p5)
rm(p6)
rm(p7)
rm(p8)
rm(image.list)

############################################################# ############################################################# #############################################################
############################################################# DE Analysis #############################################################
############################################################# ############################################################# #############################################################

seurat.object <- SCTransform(seurat.object, assay = "Spatial", verbose = T)
DefaultAssay(seurat.object) <- "SCT"
seurat.object.var <- seurat.object@assays[["SCT"]]@var.features

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = seurat.object.var, assay="Spatial")

Idents(seurat.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-0.2lnFC.txt", sep="\t")

Idents(seurat.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "integrated_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("integrated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("integrated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("integrated_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("integrated_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("integrated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)

image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")


image.list <- c("S01")
p2 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p2
ggsave("top.markers_SpatialFeaturePlots-cropped-S01.pdf", plot = p2, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S03.1")
p3 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p3
ggsave("top.markers_SpatialFeaturePlots-cropped-S03.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S04.2")
p4 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p4
ggsave("top.markers_SpatialFeaturePlots-cropped-S04.pdf", plot = p4, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S15.3")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S15.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S16.4")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S16.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S17.5")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S17.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S18.6")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S18.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("19.7")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S19.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S20.8")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S20.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S21.9")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S21.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S22.10")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S22.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S23a.11")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S23a.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S23b.12")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S23b.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S24.13")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S24.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S25.14")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S25.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S26.15")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S26.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")



top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
image.list <- c("S01")
p2 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p2
ggsave("top.markers_SpatialFeaturePlots-uncropped-S01.pdf", plot = p2, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S03.1")
p3 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p3
ggsave("top.markers_SpatialFeaturePlots-uncropped-S03.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S04.2")
p4 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p4
ggsave("top.markers_SpatialFeaturePlots-uncropped-S04.pdf", plot = p4, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S15.3")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p5
ggsave("top.markers_SpatialFeaturePlots-uncropped-S15.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S16.4")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p5
ggsave("top.markers_SpatialFeaturePlots-uncropped-S16.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S17.5")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p5
ggsave("top.markers_SpatialFeaturePlots-uncropped-S17.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S18.6")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p5
ggsave("top.markers_SpatialFeaturePlots-uncropped-S18.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("19.7")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p5
ggsave("top.markers_SpatialFeaturePlots-uncropped-S19.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S20.8")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p5
ggsave("top.markers_SpatialFeaturePlots-uncropped-S20.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S21.9")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p5
ggsave("top.markers_SpatialFeaturePlots-uncropped-S21.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S22.10")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p5
ggsave("top.markers_SpatialFeaturePlots-uncropped-S22.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S23a.11")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p5
ggsave("top.markers_SpatialFeaturePlots-uncropped-S23a.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S23b.12")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p5
ggsave("top.markers_SpatialFeaturePlots-uncropped-S23b.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S24.13")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p5
ggsave("top.markers_SpatialFeaturePlots-uncropped-S24.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S25.14")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p5
ggsave("top.markers_SpatialFeaturePlots-uncropped-S25.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S26.15")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5, crop = FALSE)
p5
ggsave("top.markers_SpatialFeaturePlots-uncropped-S26.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")




rm(de_markers)
rm(top2)
rm(top5)
rm(top5.heatmap)
rm(feature.plot)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(unique.top2)
rm(top.features)


## subset cluster 12 to identify immune cells
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(seurat.object, features = lymphoid.lineage)
png(file=paste0("lymphoid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Lymphoid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

## canonical markers for myeloid cells https://docs.abcam.com/pdf/immunology/myeloid_cell.pdf
# myeloid.lineage <- c("CD34", "CD117", "CD271", "CD33", "CD71", "CD61", "CD123", "CD44", "CD15", "CD16", "CD11b", "CD14", "CD1c", "CD141")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
myeloid.lineage <- c("CD34", "KIT", "NGFR", "CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
feature.plot <- DotPlot(seurat.object, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(seurat.object, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

macrophage.lineage <- c("CD14", "CD68", "CSF1R", "ADGRE1", "SIGLEC1", "CD38", "GPR18", "FPR2", "MYC", "EGR2")
feature.plot <- DotPlot(seurat.object, features = macrophage.lineage)
png(file=paste0("macrophage.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Macrophage Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

# https://rep.bioscientifica.com/view/journals/rep/152/5/447.xml M1, M2a, M2b, M2c
m2.lineage <- c("CD80", "CD11B", "CD36", "CD163", "CD206", "HLA-DR", "CD209", "IL1RN", "CD86", "TNFA", "IL6", "CD14")
    # following requested variables were not found: CD11B, CD206, HLA-DR, TNFA
feature.plot <- DotPlot(seurat.object, features = m2.lineage)
png(file=paste0("m2.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Macrophage Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

############################################################# ############################################################# #############################################################
############################################################# User Annotation of UMAP Clusters #############################################################
############################################################# ############################################################# #############################################################
dir.create("annotated")
setwd("annotated")
## For cluster Subtypes, I took the [integrated_DEGs_byclusters_pos-0.693lnFC.txt] 
    # cluster-Subtypes_DEGs_byclusters_top5k-pos-0.693lnFC.xlsx
  # with adj. p<0.05 and log2FC>2 and upload top marker genes into https://placentacellenrich.gdcb.iastate.edu and annotate based on vento-tormo et al or suryawanshi datasets. 
  ## Also examining the human protein atlas (https://www.proteinatlas.org), Vento-Tormo and Suryawanshi datasets (https://placentacellenrich.gdcb.iastate.edu), and the PangloaDB (https://panglaodb.se/search.html)
Idents(seurat.object) <- "seurat_clusters"
cluster.Subtypes<- c("SYT_1", "Fibro_1", "VCT_1", "Mac_1", "Pericytes", "Fibro_2", "EVT_1", "RBC", "EVT_2", "Fibro_3", "Fibro_4", "Mac_2", "SYT_2", "SYT_3", "Mac_3", "Stromal_1", "Mac_4", "Mac_5", "Mac_6", "Mac_7", "Neutrophil", "Stromal_2")
names(cluster.Subtypes) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Subtypes)
## add cluster Subtypes as a metadata column
seurat.object$Subtype <- Idents(seurat.object)

Idents(seurat.object) <- "seurat_clusters"
cluster.Type <- c("SYT", "Fibroblast", "VCT", "Macrophage", "Pericytes", "Fibroblast", "EVT", "Erythrocytes", "EVT", "Fibroblast", "Fibroblast", "Macrophage", "SYT", "SYT", "Macrophage", "Stromal", "Macrophage", "Macrophage", "Macrophage", "Macrophage", "Neutrophil", "Stromal")
names(cluster.Type) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Type)
## add cluster Types as a metadata column
seurat.object$Type <- Idents(seurat.object)


seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

## Re-order the clusters to group similar clusters together
Idents(seurat.object) <- "Subtype"
levels(seurat.object)
seurat.object2 <- subset(x = seurat.object, idents = c("SYT_1", "Fibro_1", "VCT_1", "Mac_1", "Pericytes", "Fibro_2", "EVT_1", "RBC", "EVT_2", "Fibro_3", "Fibro_4", "Mac_2", "SYT_2", "SYT_3", "Mac_3", "Stromal_1", "Mac_4", "Mac_5", "Mac_6", "Mac_7", "Neutrophil", "Stromal_2"))
ncol(seurat.object2)
    # make sure you didn't lose any cells
        # 3908 cells
levels(seurat.object2)
levels(x=seurat.object2) <- c("VCT_1", "EVT_1", "EVT_2", "SYT_1", "SYT_2", "SYT_3", "Fibro_1", "Fibro_2", "Fibro_3", "Fibro_4", "Pericytes", "Stromal_1", "Stromal_2", "RBC", "Mac_1", "Mac_2", "Mac_3", "Mac_4", "Mac_5", "Mac_6", "Mac_7", "Neutrophil")
levels(seurat.object2)
ncol(seurat.object)
  
seurat.object <- seurat.object2
rm(seurat.object2)
levels(seurat.object)

## Re-order the clusters to group similar clusters together
Idents(seurat.object) <- "Type"
levels(seurat.object)
seurat.object2 <- subset(x = seurat.object, idents = c("VCT", "EVT", "SYT", "Fibroblast", "Pericytes", "Stromal", "Erythrocytes", "Macrophage", "Neutrophil"))
ncol(seurat.object2)
        # make sure you didn't lose any cells
        # 3621 cells
levels(seurat.object2)
levels(x=seurat.object2) <- c("VCT", "EVT", "SYT", "Fibroblast", "Pericytes", "Stromal", "Erythrocytes", "Macrophage", "Neutrophil")
levels(seurat.object2)
ncol(seurat.object)
    # 3621
seurat.object <- seurat.object2
rm(seurat.object2)
levels(seurat.object)

seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
save.image("human_spatial_data-integrated-annotated_v1.RData")


Idents(seurat.object) <- "seurat_clusters"
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE")
p1
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE)
p2
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Subtype")
p3
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Condition")
umap.combined <- CombinePlots(plots = list(p4, p1, p2, p3))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("integrated_UMAP-annotated.pdf", plot = umap.combined, device = "pdf", width = 12, height = 7, units = "in")



Idents(seurat.object) <- "Type"
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
plot2<- SpatialDimPlot(seurat.object, images=image.list) vlnplot_SRY_orig.ident.pdf
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_cropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialDimPlot(seurat.object, images=image.list, crop=FALSE) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")


list.1 <- c("COL1A1", "COL3A1", "HLA-G", "IGFBP1", "PRL", "LUM")
Idents(seurat.object) <- "Type"
image.list <- c("S01", "S03.1", "S04.2")
plot2<- SpatialFeaturePlot(seurat.object, images=image.list, features=c("COL1A1", "COL3A1"), min.cutoff=0, max.cutoff=9) + patchwork::plot_layout(ncol = 3)
ggsave("Spatial-S01.markers-normal.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialFeaturePlot(seurat.object, images=image.list, features=c("HLA-G", "IGFBP1"), min.cutoff=0, max.cutoff=9) + patchwork::plot_layout(ncol = 3)
ggsave("Spatial-S03.markers-normal.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialFeaturePlot(seurat.object, images=image.list, features=c("PRL", "LUM"), min.cutoff=0, max.cutoff=9) + patchwork::plot_layout(ncol = 3)
ggsave("Spatial-S04.markers-normal.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
plot2<- SpatialFeaturePlot(seurat.object, images=image.list, features=c("ACE2"), min.cutoff=0, max.cutoff=9) + patchwork::plot_layout(ncol = 4)
ggsave("Spatial-SARS-CoV-2-receptors_ACE2", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialFeaturePlot(seurat.object, images=image.list, features=c("TMPRSS2"), min.cutoff=0, max.cutoff=9) + patchwork::plot_layout(ncol = 4)
ggsave("Spatial-SARS-CoV-2-receptors_TMPRSS2", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

Idents(seurat.object) <- "Subtype"
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
plot2<- SpatialDimPlot(seurat.object, images=image.list) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_cropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialDimPlot(seurat.object, images=image.list, crop=FALSE) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

##### Perform differential expression between Types using the raw data
subtypes <- seurat.object <- SCTransform(seurat.object, assay = "Spatial", verbose = T)
DefaultAssay(seurat.object) <- "SCT"
seurat.object.var <- seurat.object@assays[["SCT"]]@var.features
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = seurat.object.var, assay="Spatial")

Idents(seurat.object) <- "Type"
de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "Type_annotated_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("Type_annotated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("Type_annotated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Type_annotated_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Type_annotated_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Type_annotated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

##### Perform differential expression between Subtypes using the raw data
subtypes <- seurat.object <- SCTransform(seurat.object, assay = "Spatial", verbose = T)
DefaultAssay(seurat.object) <- "SCT"
seurat.object.var <- seurat.object@assays[["SCT"]]@var.features
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = seurat.object.var, assay="Spatial")

Idents(seurat.object) <- "Subtype"
de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "Subtype_annotated_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("Subtype_annotated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("Subtype_annotated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Subtype_annotated_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Subtype_annotated_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Subtype_annotated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()




setwd("..")



############################################################# ############################################################# #############################################################
############################################################# Gene set analysis #############################################################
############################################################# ############################################################# #############################################################

seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(seurat.object, features = lymphoid.lineage)
png(file=paste0("lymphoid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Lymphoid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

## canonical markers for myeloid cells https://docs.abcam.com/pdf/immunology/myeloid_cell.pdf
# myeloid.lineage <- c("CD34", "CD117", "CD271", "CD33", "CD71", "CD61", "CD123", "CD44", "CD15", "CD16", "CD11b", "CD14", "CD1c", "CD141")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
myeloid.lineage <- c("CD34", "KIT", "NGFR", "CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
feature.plot <- DotPlot(seurat.object, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(seurat.object, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

feature.plot <- DotPlot(seurat.object, features = SARS.genes)
png(file=paste0("SARS.genes-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_counts.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_logfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
rm(library.averages.heatmap)
############################################################# ############################################################# #############################################################
############################################################# DE Analysis Btwn orig.ident and uninf v inf. #############################################################
############################################################# ############################################################# #############################################################

dir.create("library")
setwd("library")
Idents(seurat.object)
Idents(seurat.object) <- "orig.ident"

p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident")
ggsave("integrated_UMAP-orig.ident.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes, assay="Spatial")

Idents(seurat.object) <- "orig.ident"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "integrated_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("integrated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("integrated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("integrated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="orig.ident")
dev.off()

image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
image.list <- c("S01")
p2 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p2
ggsave("top.markers_SpatialFeaturePlots-cropped-S01.pdf", plot = p2, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S03.1")
p3 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p3
ggsave("top.markers_SpatialFeaturePlots-cropped-S03.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S04.2")
p4 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p4
ggsave("top.markers_SpatialFeaturePlots-cropped-S04.pdf", plot = p4, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S15.3")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S15.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S16.4")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S16.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S17.5")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S17.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S18.6")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S18.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("19.7")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S19.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S20.8")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S20.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S21.9")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S21.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S22.10")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S22.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S23a.11")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S23a.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S23b.12")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S23b.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S24.13")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S24.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S25.14")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S25.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S26.15")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S26.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")



lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(seurat.object, features = lymphoid.lineage)
png(file=paste0("lymphoid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Lymphoid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

## canonical markers for myeloid cells https://docs.abcam.com/pdf/immunology/myeloid_cell.pdf
# myeloid.lineage <- c("CD34", "CD117", "CD271", "CD33", "CD71", "CD61", "CD123", "CD44", "CD15", "CD16", "CD11b", "CD14", "CD1c", "CD141")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
myeloid.lineage <- c("CD34", "KIT", "NGFR", "CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
feature.plot <- DotPlot(seurat.object, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(seurat.object, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

#seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
#save.image("human_spatial_data-integrated-annotated_v1.RData")

setwd("..")
############################################################# ############################################################# #############################################################
############################################################# Viral Transcriptomics #############################################################
############################################################# ############################################################# #############################################################
dir.create("viral.transcriptomics")
setwd("viral.transcriptomics")
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")


## kernel density estimation 
## https://www.bioconductor.org/packages/release/bioc/vignettes/Nebulosa/inst/doc/nebulosa_seurat.html#1_Overview
  ## https://academic.oup.com/bioinformatics/article-abstract/37/16/2485/6103785
# BiocManager::install("Nebulosa")
library("Nebulosa")
# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
p1 <- plot_density(seurat.object, "ORF1AB-1")
p2 <- FeaturePlot(seurat.object, "ORF1AB-1")
p5 <- CombinePlots(plots = list(p1, p2))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("integrated_KernelDensity_ORF1AB-1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots
p3 <- plot_density(seurat.object, c("ORF1AB-1", "ORF10"))
p5,- p3 + plot_layout(ncol = 1)
ggsave("integrated_KernelDensity_Joint_ORF1AB-1_ORF10.pdf", plot = p5, device = "pdf", width = 5, height = 7, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("ORF1AB-1", "ORF10"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("integrated_KernelDensity_Joint_ORF1AB-1_ORF10.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")



image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
plot2<- SpatialFeaturePlot(seurat.object, images=image.list, features=c("ACE2"), min.cutoff=0, max.cutoff=9) + patchwork::plot_layout(ncol = 4)
ggsave("Spatial-SARS-CoV-2-receptors_ACE2", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialFeaturePlot(seurat.object, images=image.list, features=c("TMPRSS2"), min.cutoff=0, max.cutoff=9) + patchwork::plot_layout(ncol = 4)
ggsave("Spatial-SARS-CoV-2-receptors_TMPRSS2", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")


ncol(seurat.object) # 17927 total spots
acute.viral <- subset(seurat.object, percent.viral > 0)

ax.acute.viral <- subset(seurat.object, ident= "SARS-CoV-2.ax")
ncol(ax.acute.viral)  # 5244
sx.acute.viral <- subset(seurat.object, ident= "SARS-CoV-2.sx")
ncol(sx.acute.viral)  # 4202
uninfected.acute.viral <- subset(seurat.object, ident= "Uninfected")
ncol(uninfected.acute.viral)  # 8481 

> levels(acute.viral)
# [1] "SARS-CoV-2.sx"
ax.acute.viral <- subset(acute.viral, ident= "SARS-CoV-2.ax")
ncol(ax.acute.viral)  # 5244
sx.acute.viral <- subset(acute.viral, ident= "SARS-CoV-2.sx")
ncol(sx.acute.viral)  # 752 out of 4202 spots with viral transcripts
uninfected.acute.viral <- subset(acute.viral, ident= "Uninfected")
ncol(uninfected.acute.viral)  # 8481 ??

rm(uninfected.acute.viral)
rm(ax.acute.viral)
rm(sx.acute.viral)

ncol(acute.viral)   # 752 with viral transcripts
acute.viral <- subset(acute.viral, features = rownames(acute.viral)[rownames(acute.viral) %in% SARS.genes])
write.table(acute.viral@assays[["Spatial"]]@counts, "viral.counts_infected-cells.txt", sep="\t")
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 4)
ggsave("infected-only_SpatialFeaturePlot_percent.viral.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

# image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
image.list <- c("S01")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S01_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S03.1")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S03_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S04.2")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S04_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S15.3")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S15_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S16.4")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S16_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S17.5")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S17_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S18.6")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S18_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S19.7")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S19_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S20.8")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S20_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S21.9")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S21_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S22.10")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S22_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23a.11")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23a_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23b.12")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23b_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S24.13")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S24_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S25.14")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S25_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S26.15")
plot2<- SpatialFeaturePlot(seurat.object, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S26_SpatialFeaturePlot_percent.viral.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")








seurat.object <- ScaleData(seurat.object, features = SARS.genes, assay = "Spatial")


Idents(seurat.object) <- "Type"
feature.plot <- DotPlot(seurat.object, features = SARS.genes)
png(file=paste0("SARS.genes-DotPlot_Type.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_counts_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_logfc_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(seurat.object) <- "Condition"
feature.plot <- DotPlot(seurat.object, features = SARS.genes)
png(file=paste0("SARS.genes-DotPlot_Condition.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="Condition")
dev.off()
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_counts_Condition.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_logfc_Condition.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(seurat.object) <- "orig.ident"
feature.plot <- DotPlot(seurat.object, features = SARS.genes)
png(file=paste0("SARS.genes-DotPlot_Condition.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="orig.ident")
dev.off()
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")


SARS.receptors <- c("ACE2", "TMPRSS2", "CTSB", "CTSL")
SARS.genes.receptors <- union(SARS.genes, SARS.receptors)
seurat.object <- ScaleData(seurat.object, features = SARS.genes.receptors, assay = "Spatial")

library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes.receptors, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors_seurat.clusters_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes.receptors, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors_seurat.clusters_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes.receptors, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors_infected.only_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes.receptors, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors_infected.only_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

library.averages.heatmap <- VlnPlot(acute.viral, features = "ACE2", assay="Spatial")
ggsave("SARS.genes.receptors_infected.only_vlnplot_ACE2_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- VlnPlot(acute.viral, features = "TMPRSS2", assay="Spatial")
ggsave("SARS.genes.receptors_infected.only_vlnplot_TMPRSS2_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(acute.viral) <- "orig.ident"
library.averages.heatmap <- VlnPlot(acute.viral, features = "SRY", assay="Spatial")
ggsave("SARS.genes.receptors_infected.only_vlnplot_SRY_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

sex.genes <- c("SRY", "DDX3Y")

Idents(seurat.object) <- "orig.ident"
library.averages.heatmap <- VlnPlot(seurat.object, features = "SRY", assay="Spatial")
ggsave("vlnplot_SRY_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object) <- "orig.ident"
library.averages.heatmap <- VlnPlot(seurat.object, features = "XIST", assay="Spatial")
ggsave("vlnplot_XIST_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object) <- "orig.ident"
library.averages.heatmap <- VlnPlot(seurat.object, features = "HBB", assay="Spatial")
ggsave("vlnplot_HBB_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object) <- "orig.ident"
library.averages.heatmap <- VlnPlot(seurat.object, features = "DDX3X", assay="Spatial")
ggsave("vlnplot_DDX3X_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object) <- "orig.ident"
library.averages.heatmap <- VlnPlot(seurat.object, features = "DDX3Y", assay="Spatial")
ggsave("vlnplot_DDX3Y_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

DDX3Y.genes <- c("DDX3Y")
seurat.object <- PercentageFeatureSet(seurat.object, features = DDX3Y.genes, col.name = "percent.DDX3Y")
seurat.object.y <- subset(seurat.object, percent.DDX3Y > 0)
ncol(seurat.object.y)   # 3499

image.list <- c("S16.4")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.DDX3Y", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S16_SpatialFeaturePlot_percent.DDX3Y.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S18.6")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.DDX3Y", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S18_SpatialFeaturePlot_percent.DDX3Y.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S19.7")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.DDX3Y", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S19_SpatialFeaturePlot_percent.DDX3Y.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S22.10")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.DDX3Y", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S22_SpatialFeaturePlot_percent.DDX3Y.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S24.13")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.DDX3Y", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S24_SpatialFeaturePlot_percent.DDX3Y.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S25.14")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.DDX3Y", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S25_SpatialFeaturePlot_percent.DDX3Y.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S26.15")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.DDX3Y", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S26_SpatialFeaturePlot_percent.DDX3Y.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")


XIST.genes <- c("XIST")
seurat.object <- PercentageFeatureSet(seurat.object, features = XIST.genes, col.name = "percent.XIST")
seurat.object.y <- subset(seurat.object, percent.XIST < 0.00000000000000000001)
ncol(seurat.object.y)   # 12503 12503

image.list <- c("S01")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S01_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S03.1")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S03_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S04.2")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S04_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S15.3")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S15_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S16.4")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S16_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S17.5")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S17_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S18.6")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S18_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S19.7")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S19_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S20.8")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S20_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S21.9")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S21_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S22.10")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S22_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23a.11")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23a_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23b.12")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23b_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S24.13")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S24_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S25.14")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S25_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S26.15")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.XIST", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S26_SpatialFeaturePlot_percent.XIST.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")


ACE2.genes <- c("ACE2")
seurat.object <- PercentageFeatureSet(seurat.object, features = ACE2.genes, col.name = "percent.ACE2")
seurat.object.y <- subset(seurat.object, percent.ACE2 > 0)
ncol(seurat.object.y)   # 391
acute.viral <- subset(seurat.object.y, percent.viral > 0)
ncol(acute.viral)   # 0

image.list <- c("S01")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S01_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S03.1")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S03_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S04.2")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S04_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S15.3")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S15_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S16.4")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S16_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S17.5")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S17_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S18.6")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S18_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S19.7")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S19_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S20.8")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S20_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S21.9")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S21_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S22.10")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S22_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23a.11")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23a_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23b.12")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23b_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S24.13")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S24_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S25.14")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S25_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S26.15")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S26_SpatialFeaturePlot_percent.ACE2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")

image.list <- c("S01")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S01_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S03.1")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S03_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S04.2")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S04_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S15.3")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S15_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S16.4")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S16_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S17.5")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S17_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S18.6")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S18_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S19.7")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S19_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S20.8")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S20_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S21.9")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S21_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S22.10")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S22_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23a.11")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23a_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23b.12")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23b_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S24.13")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S24_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S25.14")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S25_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S26.15")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S26_SpatialFeaturePlot_percent.ACE2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")



TMPRSS2.genes <- c("TMPRSS2")
seurat.object <- PercentageFeatureSet(seurat.object, features = TMPRSS2.genes, col.name = "percent.TMPRSS2")
seurat.object.y <- subset(seurat.object, percent.TMPRSS2 > 0)
ncol(seurat.object.y)   # 59
acute.viral <- subset(seurat.object.y, percent.viral > 0)
ncol(acute.viral)   # 0

image.list <- c("S01")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S01_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S03.1")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S03_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S04.2")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S04_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S15.3")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S15_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S16.4")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S16_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S17.5")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S17_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S18.6")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S18_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S19.7")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S19_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S20.8")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S20_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S21.9")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S21_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S22.10")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S22_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23a.11")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23a_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23b.12")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23b_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S24.13")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S24_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S25.14")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S25_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S26.15")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S26_SpatialFeaturePlot_percent.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")


image.list <- c("S01")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S01_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S03.1")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S03_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S04.2")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S04_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S15.3")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S15_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S16.4")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S16_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S17.5")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S17_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S18.6")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S18_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S19.7")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S19_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S20.8")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S20_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S21.9")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S21_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S22.10")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S22_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23a.11")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23a_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23b.12")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23b_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S24.13")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S24_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S25.14")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S25_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S26.15")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S26_SpatialFeaturePlot_percent.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")

ACE2.TMPRSS2.genes <- c("ACE2", "TMPRSS2")
seurat.object <- PercentageFeatureSet(seurat.object, features = ACE2.TMPRSS2.genes, col.name = "percent.ACE2.TMPRSS2")
seurat.object.y <- subset(seurat.object, percent.ACE2.TMPRSS2 > 0)
ncol(seurat.object.y)   # 446
acute.viral <- subset(seurat.object.y, percent.viral > 0)
ncol(acute.viral)   # 0
acute.viral <- subset(seurat.object.y, percent.ACE2 < 0.0000000001)
ncol(acute.viral)   # 55
acute.viral <- subset(seurat.object.y, percent.TMPRSS2 < 0.0000000001)
ncol(acute.viral)   # 387
acute.viral <- subset(acute.viral, percent.ACE2 < 0.0000000001)
    # Error: No cells found

image.list <- c("S01")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S01_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S03.1")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S03_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S04.2")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S04_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S15.3")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S15_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S16.4")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S16_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S17.5")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S17_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S18.6")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S18_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S19.7")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S19_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S20.8")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S20_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S21.9")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S21_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S22.10")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S22_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23a.11")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23a_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23b.12")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23b_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S24.13")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S24_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S25.14")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S25_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S26.15")
plot2<- SpatialFeaturePlot(seurat.object.y, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S26_SpatialFeaturePlot_percent.ACE2.TMPRSS2.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")


image.list <- c("S01")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S01_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S03.1")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S03_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S04.2")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S04_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S15.3")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S15_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S16.4")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S16_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S17.5")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S17_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S18.6")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S18_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S19.7")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S19_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S20.8")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S20_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S21.9")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S21_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S22.10")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S22_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23a.11")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23a_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S23b.12")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S23b_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S24.13")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S24_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S25.14")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S25_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")
image.list <- c("S26.15")
plot2<- SpatialFeaturePlot(acute.viral, features = "percent.ACE2.TMPRSS2", images=image.list) + patchwork::plot_layout(ncol = 1)
ggsave("S26_SpatialFeaturePlot_percent.ACE2.TMPRSS2inf-only.png", plot = plot2, device = "png", width = 8, height = 11, units = "in")













acute.viral.nobld <- subset(acute.viral, percent.HBB < 0.000000001)
ncol(acute.viral.nobld) # 710
library.averages.heatmap <- VlnPlot(acute.viral.nobld, features = "HBB", assay="Spatial")
ggsave("vlnplot_HBB_inf.only-noHBB_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

SARS.genes.receptors.HBB <- union(HBB.genes, SARS.genes.receptors)
acute.viral.nobld <- ScaleData(acute.viral.nobld, features = SARS.genes.receptors.HBB, assay = "Spatial")


library.averages.heatmap <- DoHeatmap(acute.viral.nobld, features = SARS.genes.receptors.HBB, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors.HBB_infected.only-nobld_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(acute.viral.nobld, features = SARS.genes.receptors.HBB, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors.HBB_infected.only-nobld_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(acute.viral.nobld) <- "Type"
library.averages.heatmap <- DoHeatmap(acute.viral.nobld, features = SARS.genes.receptors.HBB, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors.HBB_infected.only-nobld_heatmap_counts_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(acute.viral.nobld, features = SARS.genes.receptors.HBB, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors.HBB_infected.only-nobld_heatmap_logfc_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

acute.viral.nobld <- subset(acute.viral, percent.Heme < 0.000000001)
ncol(acute.viral.nobld) # 611

SARS.genes.receptors.HBB <- union(HB.genes, SARS.genes.receptors)
acute.viral.nobld <- ScaleData(acute.viral.nobld, features = SARS.genes.receptors.HBB, assay = "Spatial")

Idents(acute.viral.nobld) <- "orig.ident"
library.averages.heatmap <- DoHeatmap(acute.viral.nobld, features = SARS.genes.receptors.HBB, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors.HBB_infected.only-noheme_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(acute.viral.nobld, features = SARS.genes.receptors.HBB, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors.HBB_infected.only-noheme_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(acute.viral.nobld) <- "Type"
library.averages.heatmap <- DoHeatmap(acute.viral.nobld, features = SARS.genes.receptors.HBB, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors.HBB_infected.only-noheme_heatmap_counts_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(acute.viral.nobld, features = SARS.genes.receptors.HBB, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors.HBB_infected.only-noheme_heatmap_logfc_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

acute.viral <- ScaleData(acute.viral, features = SARS.genes.receptors.HBB, assay = "Spatial")
Idents(acute.viral) <- "orig.ident"
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes.receptors.HBB, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors.HBB_infected.only_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes.receptors.HBB, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors.HBB_infected.only_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(acute.viral) <- "Type"
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes.receptors.HBB, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors.HBB_infected.only_heatmap_counts_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes.receptors.HBB, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors.HBB_infected.only_heatmap_logfc_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")



acute.viral <- ScaleData(acute.viral, features = SARS.genes, assay = "Spatial")

Idents(acute.viral) <- "Type"
feature.plot <- DotPlot(acute.viral, features = SARS.genes)
png(file=paste0("SARS.genes-DotPlot_Typeinfected-only.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_counts_Typeinfected-only_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_logfc_Typeinfected-only_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 8, units = "in")
acute.viral <- ScaleData(acute.viral, features = SARS.genes.receptors, assay = "Spatial")
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes.receptors, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors_seurat.clusters_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes.receptors, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors_seurat.clusters_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")


Idents(acute.viral) <- "Subtype"
feature.plot <- DotPlot(acute.viral, features = SARS.genes)
png(file=paste0("SARS.genes-DotPlot_Subtypeinfected-only.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_counts_Subtypeinfected-only_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_logfc_Subtypeinfected-only_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 8, units = "in")




Idents(acute.viral) <- "Condition"
feature.plot <- DotPlot(acute.viral, features = SARS.genes)
png(file=paste0("SARS.genes-DotPlot_Conditioninfected-only.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="Condition")
dev.off()
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_counts_Conditioninfected-only_Condition.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_logfc_Conditioninfected-only_Condition.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 8, units = "in")

Idents(acute.viral) <- "orig.ident"
feature.plot <- DotPlot(acute.viral, features = SARS.genes)
png(file=paste0("SARS.genes-DotPlot_Conditioninfected-only.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="orig.ident")
dev.off()
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_counts_Conditioninfected-only_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(acute.viral, features = SARS.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_logfc_Conditioninfected-only_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 8, units = "in")

rm(acute.viral)

image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")


image.list <- c("S01")
p2 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5) 
ggsave("viral.markers_SpatialFeaturePlots-cropped-S01.pdf", plot = p2, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S03.1")
p3 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S03.pdf", plot = p3, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S04.2")
p4 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S04.pdf", plot = p4, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S15.3")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S15.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S16.4")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S16.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S17.5")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S17.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S18.6")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S18.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("19.7")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S19.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S20.8")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S20.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S21.9")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S21.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S22.10")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S22.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S23a.11")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S23a.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S23b.12")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S23b.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S24.13")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S24.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S25.14")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S25.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S26.15")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S26.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")

image.list <- c("S01")
p2 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S01scaleyfixey.pdf", plot = p2, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S03.1")
p3 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S03scaleyfixey.pdf", plot = p3, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S04.2")
p4 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S04scaleyfixey.pdf", plot = p4, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S15.3")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S15scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S16.4")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S16scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S17.5")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S17scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S18.6")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S18scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("19.7")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S19scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S20.8")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S20scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S21.9")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S21scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S22.10")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S22scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S23a.11")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S23ascaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S23b.12")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S23bscaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S24.13")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S24scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S25.14")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S25scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S26.15")
p5 <- SpatialFeaturePlot(seurat.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S26scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")

seurat.object <- ScaleData(seurat.object, features = SARS.genes, assay = "Spatial")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF1AB-1", images=image.list, ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF1AB-1_SpatialFeatureplot-cropped.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "S", images=image.list, ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("S_SpatialFeatureplot-cropped.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF3A", images=image.list, ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF3A_SpatialFeatureplot-cropped.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "E", images=image.list, ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("E_SpatialFeatureplot-cropped.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "M", images=image.list, ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("M_SpatialFeatureplot-cropped.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF6", images=image.list, ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF6_SpatialFeatureplot-cropped.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF7A", images=image.list, ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF7A_SpatialFeatureplot-cropped.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF7B", images=image.list, ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF7B_SpatialFeatureplot-cropped.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF8", images=image.list, ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF8_SpatialFeatureplot-cropped.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "N", images=image.list, ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("N_SpatialFeatureplot-cropped.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF10", images=image.list, ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF10_SpatialFeatureplot-cropped.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")

# SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF1AB-1", images=image.list,min.cutoff=0, max.cutoff=9,  ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF1AB-1_SpatialFeatureplot-cropped-scaley-fixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "S", images=image.list,min.cutoff=0, max.cutoff=9,  ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("S_SpatialFeatureplot-cropped-scaley-fixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF3A", images=image.list,min.cutoff=0, max.cutoff=9,  ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF3A_SpatialFeatureplot-cropped-scaley-fixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "E", images=image.list,min.cutoff=0, max.cutoff=9,  ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("E_SpatialFeatureplot-cropped-scaley-fixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "M", images=image.list,min.cutoff=0, max.cutoff=9,  ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("M_SpatialFeatureplot-cropped-scaley-fixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF6", images=image.list,min.cutoff=0, max.cutoff=9,  ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF6_SpatialFeatureplot-cropped-scaley-fixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF7A", images=image.list,min.cutoff=0, max.cutoff=9,  ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF7A_SpatialFeatureplot-cropped-scaley-fixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF7B", images=image.list,min.cutoff=0, max.cutoff=9,  ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF7B_SpatialFeatureplot-cropped-scaley-fixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF8", images=image.list,min.cutoff=0, max.cutoff=9,  ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF8_SpatialFeatureplot-cropped-scaley-fixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "N", images=image.list,min.cutoff=0, max.cutoff=9,  ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("N_SpatialFeatureplot-cropped-scaley-fixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "ORF10", images=image.list,min.cutoff=0, max.cutoff=9,  ncol=4) + patchwork::plot_layout(ncol = 4)
ggsave("ORF10_SpatialFeatureplot-cropped-scaley-fixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")


seurat.object <- ScaleData(seurat.object, features = SARS.genes, assay = "Spatial")

Idents(seurat.object) <- "Condition"
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=SARS.genes)
write.csv(cluster.averages@assays[["Spatial"]]@counts, file = "library_viralmeancounts.csv")
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("viral_mc_heatmap_Condition.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("viral_FC_heatmap_Condition.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
SARS.genes.vlnplot <- VlnPlot(seurat.object, features = "percent.viral", pt.size = 0.1) + NoLegend()
ggsave("SARS.genes_vlnplot_Condition.png", plot = SARS.genes.vlnplot, device = "png", width = 7, height = 5, units = "in")
Idents(seurat.object) <- "seurat_clusters"

Idents(seurat.object) <- "orig.ident"
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=SARS.genes)
write.csv(cluster.averages@assays[["Spatial"]]@counts, file = "library_viralmeancounts.csv")
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("viral_mc_heatmap_orig.ident.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("viral_FC_heatmap_orig.ident.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
SARS.genes.vlnplot <- VlnPlot(seurat.object, features = "percent.viral", pt.size = 0.1) + NoLegend()
ggsave("SARS.genes_vlnplot_orig.ident.png", plot = SARS.genes.vlnplot, device = "png", width = 7, height = 5, units = "in")


Idents(seurat.object) <- "Type"
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=SARS.genes)
write.csv(cluster.averages@assays[["Spatial"]]@counts, file = "library_viralmeancounts.csv")
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("viral_mc_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("viral_FC_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
SARS.genes.vlnplot <- VlnPlot(seurat.object, features = "percent.viral", pt.size = 0.1) + NoLegend()
ggsave("SARS.genes_vlnplot_Type.png", plot = SARS.genes.vlnplot, device = "png", width = 7, height = 5, units = "in")

Idents(seurat.object) <- "seurat_clusters"


image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
plot2<- SpatialFeaturePlot(seurat.object, images=image.list, features=c("ACE2"), min.cutoff=0, max.cutoff=9) + patchwork::plot_layout(ncol = 4)
ggsave("Spatial-SARS-CoV-2-receptors_ACE2", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialFeaturePlot(seurat.object, images=image.list, features=c("TMPRSS2"), min.cutoff=0, max.cutoff=9) + patchwork::plot_layout(ncol = 4)
ggsave("Spatial-SARS-CoV-2-receptors_TMPRSS2", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")


gene.list <- ("HBB", "HBA", "SRY", "XIST", "ACE2", "TMPRSS2", "ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
Idents(seurat.object) <- "orig.ident"
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=SARS.genes)
write.csv(cluster.averages@assays[["Spatial"]]@counts, file = "gene.list_meancounts.csv")
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, slot="counts", assay="Spatial")
ggsave("gene.list_mc_heatmap_orig.ident.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = gene.list, raster = FALSE, assay="Spatial")
ggsave("gene.list_FC_heatmap_orig.ident.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


SARS.receptors <- c("ACE2", "TMPRSS2", "CTSB", "CTSL", "HBB", "HBA", "SRY", "XIST","ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes.receptors <- union(SARS.genes, SARS.receptors)

Idents(seurat.object) <- "orig.ident"
seurat.object <- ScaleData(seurat.object, features = SARS.genes.receptors, assay = "Spatial")
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes.receptors, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors-HB-Sex_orig.ident_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes.receptors, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors-HB-Sex_orig.ident_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(seurat.object) <- "Condition"
seurat.object <- ScaleData(seurat.object, features = SARS.genes.receptors, assay = "Spatial")
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes.receptors, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors-HB-Sex_Condition_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes.receptors, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors-HB-Sex_Condition_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

HB.genes <- c("HBA1", "HBB", "HBD", "HBZ", "HBQ1", "ICAM4")



HBB.genes <- c("HBB")
SRY.genes <- c("SRY")
XIST.genes <- c("XIST")
ACE2.genes <- c("ACE2")
TMPRSS2.genes <- c("TMPRSS2")
ACE2.TMPRSS2.genes <- c("ACE2", "TMPRSS2")

acute.viral <- PercentageFeatureSet(acute.viral, features = HBB.genes, col.name = "percent.HBB")
acute.viral <- PercentageFeatureSet(acute.viral, features = HB.genes, col.name = "percent.Heme")


seurat.object <- PercentageFeatureSet(seurat.object, features = HBB.genes, col.name = "percent.HBB")
seurat.object <- PercentageFeatureSet(seurat.object, features = SRY.genes, col.name = "percent.SRY")
seurat.object <- PercentageFeatureSet(seurat.object, features = XIST.genes, col.name = "percent.XIST")
seurat.object <- PercentageFeatureSet(seurat.object, features = ACE2.genes, col.name = "percent.ACE2")
seurat.object <- PercentageFeatureSet(seurat.object, features = TMPRSS2.genes, col.name = "percent.TMPRSS2")
seurat.object <- PercentageFeatureSet(seurat.object, features = TMPRSS2.genes, col.name = "percent.ACE2.TMPRSS2")


seurat.object2 <- subset(seurat.object, subset = percent.SRY > 0 )
ncol(seurat.object2)
levels(seurat.object2)
# new.metadata <- c("M", "M","M","M","M","M",)
# names(new.metadata) <- levels(seurat.object2)
# seurat.object2 <- RenameIdents(seurat.object2, new.metadata)
seurat.object2$FetalSex <-  RenameIdents(seurat.object2, new.metadata)

seurat.object3 <- subset(seurat.object, subset = percent.SRY < 0 )
seurat.object4 <- merge(seurat.object3, seurat.object4)


seurat.object2 <- subset(seurat.object, subset = percent.ACE2 > 0 )
ncol(seurat.object2)
seurat.object2 <- subset(seurat.object2, subset = percent.viral > 0 )
ncol(seurat.object2)
SARS.receptors <- c("ACE2", "TMPRSS2", "CTSB", "CTSL", "HBB", "HBA", "SRY", "XIST","ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes.receptors <- union(SARS.genes, SARS.receptors)
Idents(seurat.object) <- "orig.ident"
seurat.object <- ScaleData(seurat.object2, features = SARS.genes.receptors, assay = "Spatial")
library.averages.heatmap <- DoHeatmap(seurat.object2, features = SARS.genes.receptors, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("Ace2.SARS.genes.receptors-HB-Sex_orig.ident_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object2, features = SARS.genes.receptors, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("Ace2.SARS.genes.receptors-HB-Sex_orig.ident_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")


seurat.object2 <- subset(seurat.object, subset = percent.TMPRSS2 > 0 )
ncol(seurat.object2)
seurat.object2 <- subset(seurat.object2, subset = percent.viral > 0 )
ncol(seurat.object2)
SARS.genes.receptors <- union(SARS.genes, SARS.receptors)
Idents(seurat.object) <- "orig.ident"
seurat.object <- ScaleData(seurat.object2, features = SARS.genes.receptors, assay = "Spatial")
library.averages.heatmap <- DoHeatmap(seurat.object2, features = SARS.genes.receptors, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("Tmprss2.SARS.genes.receptors-HB-Sex_orig.ident_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object2, features = SARS.genes.receptors, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("Tmprss2.SARS.genes.receptors-HB-Sex_orig.ident_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")


seurat.object2 <- subset(seurat.object, subset = percent.ACE2.TMPRSS2 > 0 )
ncol(seurat.object2)
seurat.object2 <- subset(seurat.object2, subset = percent.viral > 0 )
ncol(seurat.object2)
SARS.genes.receptors <- union(SARS.genes, SARS.receptors)
Idents(seurat.object) <- "orig.ident"
seurat.object <- ScaleData(seurat.object2, features = SARS.genes.receptors, assay = "Spatial")
library.averages.heatmap <- DoHeatmap(seurat.object2, features = SARS.genes.receptors, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("ace2.Tmprss2.SARS.genes.receptors-HB-Sex_orig.ident_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object2, features = SARS.genes.receptors, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("ace2.Tmprss2.SARS.genes.receptors-HB-Sex_orig.ident_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")


levels(seurat.object2)
# new.metadata <- c("M", "M","M","M","M","M",)
# names(new.metadata) <- levels(seurat.object2)
# seurat.object2 <- RenameIdents(seurat.object2, new.metadata)
seurat.object2$FetalSex <-  RenameIdents(seurat.object2, new.metadata)



setwd("..")
############################################################# ############################################################# #############################################################
############################################################# ASSESSMENT OF KNOWN CONDTIONS (e.g. uninfected v infected) #############################################################
############################################################# ############################################################# #############################################################
############ Condition

dir.create("Condition")
setwd("Condition")
# Idents(seurat.object) <- "orig.ident"
# levels(seurat.object)
    # "S01"  "S03"  "S04"  "S15"  "S16"  "S17"  "S18"  "S19"  "S20"  "S21"  "S22"  "S23a" "S23b" "S24"  "S25"  "S26" 
  # Metadata
        ## uninfected: S01, S02, S03, S24, S25, S26
        ## SARS-CoV-2 asymptomatic: S15, S17, S18, S19
        ## SARS-CoV-2 symptomatic: S16, S20, S21, S22, S23a, S23b
          ## SARS-CoV-2 perivillous fibrinoids: S23a and S23b
          ## Noted areas of inflammation: S18
        ## Samples S01, S02, and S03 were from the sample placenta from distinct regions including chorionic villi, decidua, and a chorioamiotic membrane roll
        ## Samples 15-26 were sampled from the placental parenchyme with various regions of villi and some decidua
        ## Samples 23a and 23b were from the same placenta punch, selected for two sections due to the presence of perivillious fibrinoid pathology identified in QC H&E slide


# new.metadata <- c("Uninfected", "Uninfected", "Uninfected", "SARS-CoV-2.ax", "SARS-CoV-2.sx", "SARS-CoV-2.ax", "SARS-CoV-2.ax", "SARS-CoV-2.ax", "SARS-CoV-2.sx", "SARS-CoV-2.sx", "SARS-CoV-2.sx", "SARS-CoV-2.sx", "SARS-CoV-2.sx", "Uninfected","Uninfected","Uninfected")
# names(new.metadata) <- levels(seurat.object)
# seurat.object <- RenameIdents(seurat.object, new.metadata)

## add cluster Subtypes as a metadata column
# seurat.object$Condition <- Idents(seurat.object)

Idents(seurat.object) <- "Condition"
Idents(seurat.object)
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident")
ggsave("integrated_UMAP-Condition.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "Condition_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "Condition_DEGs_byclusters_pos-0.2lnFC.txt", sep="\t")

seurat.object <- ScaleData(seurat.object, features = all.genes)

top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("Condition_top15_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("Condition_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Condition_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Condition")
dev.off()
rm(new.metadata)
rm(feature.plot)
rm(top5.heatmap)

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)

image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
image.list <- c("S01")
p2 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p2
ggsave("top.markers_SpatialFeaturePlots-cropped-S01.pdf", plot = p2, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S03.1")
p3 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p3
ggsave("top.markers_SpatialFeaturePlots-cropped-S03.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S04.2")
p4 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p4
ggsave("top.markers_SpatialFeaturePlots-cropped-S04.pdf", plot = p4, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S15.3")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S15.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S16.4")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S16.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S17.5")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S17.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S18.6")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S18.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("19.7")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S19.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S20.8")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S20.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S21.9")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S21.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S22.10")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S22.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S23a.11")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S23a.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S23b.12")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S23b.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S24.13")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S24.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S25.14")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S25.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S26.15")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S26.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")


setwd("..")


############################################################# ############################################################# #############################################################
############################################################# subset macrophages from human spatial data #############################################################
############################################################# ############################################################# #############################################################

set.seed(seed=1)
setwd("/home/ebarrozo/visium/results")
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
library(SeuratDisk)
library(dplyr)
library(patchwork)

setwd("/home/ebarrozo/visium/results/seurat_human_v2")
load("/home/ebarrozo/visium/results/seurat_human_v2/annotated/human_spatial_data-integrated-annotated_v1.RData")
setwd("/home/ebarrozo/visium/results/seurat_human_v2")

dir.create("subset.macs")
setwd("subset.macs")
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE)


setwd("/home/ebarrozo/visium/results/seurat_human_v2/subset.macs")
Idents(seurat.object) <- "Type"
levels(seurat.object)
ncol(seurat.object)
  # 17927

seurat.object2 <- subset(x = seurat.object, idents = c("Macrophage"))
ncol(seurat.object2)
  # 3180
seurat.object <- seurat.object2
rm(seurat.object2)

############################################################# ############################################################# #############################################################
############################################################# Integrated: Dimension Reduction and Visualization #############################################################
############################################################# ############################################################# #############################################################

seurat.object <- ScaleData(seurat.object, features = integrated.genes, assay = "integrated", vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"))
seurat.object<- RunPCA(seurat.object, features = integrated.genes, assay = "integrated", slot = "scale.data")
elbow.plot <- ElbowPlot(seurat.object)
ggsave("seurat.object_elbow_plot.pdf", plot = elbow.plot, device = "pdf")
rm(elbow.plot)

DefaultAssay(seurat.object) <- "integrated"
seurat.object <- FindNeighbors(seurat.object, features = "integrated.genes", dims = 1:30)
seurat.object <- FindClusters(seurat.object, resolution = 1.2)
seurat.object <- RunUMAP(seurat.object, dims=1:30)
  ## clusters 0 thru 8

Idents(seurat.object) <- "seurat_clusters"
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE")
p1
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Condition", label= "TRUE")
p2
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Phase", label= "TRUE")
p3
p4 <- FeaturePlot(seurat.object, reduction = "umap", features = "percent.viral")
p4

umap.combined <- CombinePlots(plots = list(p1, p2, p3, p4))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("integrated_UMAP-combined.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("integrated_UMAP-combined-wide.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "orig.ident", label=TRUE)
ggsave("integrated_UMAP_splitby_libID.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", split.by = "orig.ident", label=TRUE)
ggsave("integrated_UMAP_splitby_libID_Type.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Subtype", split.by = "orig.ident", label=TRUE)
ggsave("integrated_UMAP_splitby_libID_Subtype.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "seurat_clusters", label=TRUE)
ggsave("integrated_UMAP_splitby_clusters.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
#examine UMAPs with qc metrics 
p5 <- FeaturePlot(seurat.object, features = 'nFeature_Spatial')
p6 <- FeaturePlot(seurat.object, features = 'nCount_Spatial')
p7 <- FeaturePlot(seurat.object, features = 'percent.mt')
p8 <- FeaturePlot(seurat.object, features = 'percent.ribo')

umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8), ncol=2)
ggsave("integrated_UMAP_QCmetricsFeaturePlots.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")


Idents(seurat.object) <- "seurat_clusters"
## subset cluster 12 to identify immune cells
DefaultAssay(seurat.object) <- "SCT"
seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "SCT")

## use spatial for all DE analysis and plots
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")


## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(seurat.object, features = lymphoid.lineage)
png(file=paste0("lymphoid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Lymphoid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

## canonical markers for myeloid cells https://docs.abcam.com/pdf/immunology/myeloid_cell.pdf
# myeloid.lineage <- c("CD34", "CD117", "CD271", "CD33", "CD71", "CD61", "CD123", "CD44", "CD15", "CD16", "CD11b", "CD14", "CD1c", "CD141")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
myeloid.lineage <- c("CD34", "KIT", "NGFR", "CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
feature.plot <- DotPlot(seurat.object, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(seurat.object, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

macrophage.lineage <- c("CD14", "CD68", "CSF1R", "ADGRE1", "SIGLEC1", "CD38", "GPR18", "FPR2", "MYC", "EGR2")
feature.plot <- DotPlot(seurat.object, features = macrophage.lineage)
png(file=paste0("macrophage.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Macrophage Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()


# 65. Aagaard-Tillery KM, Silver R, Dalton J. Immunology of normal pregnancy. Semin Fetal Neonatal Med. Oct 2006;11(5):279-95. doi:10.1016/j.siny.2006.04.003
# 66. Ivashkiv LB. Epigenetic regulation of macrophage polarization and function. Trends Immunol. May 2013;34(5):216-23. doi:10.1016/j.it.2012.11.001
# 67. Murray PJ. Macrophage Polarization. Annu Rev Physiol. 02 10 2017;79:541-566. doi:10.1146/annurev-physiol-022516-034339
# 68. Yao Y, Xu XH, Jin L. Macrophage Polarization in Physiological and Pathological Pregnancy. Front Immunol. 2019;10:792. doi:10.3389/fimmu.2019.00792
# 69. Mues B, Langer D, Zwadlo G, Sorg C. Phenotypic characterization of macrophages in human term placenta. Immunology. Jul 1989;67(3):303-7. 
# 70. Bulmer JN, Johnson PM. Macrophage populations in the human placenta and amniochorion. Clin Exp Immunol. Aug 1984;57(2):393-403. 
# 71. Loegl J, Hiden U, Nussbaumer E, et al. Hofbauer cells of M2a, M2b and M2c polarization may regulate feto-placental angiogenesis. Reproduction. 2016;152(5):447-455. doi:10.1530/REP-16-0159
# 74. Schliefsteiner C, Ibesich S, Wadsack C. Placental Hofbauer Cell Polarization Resists Inflammatory Cues In Vitro. Int J Mol Sci. Jan 22 2020;21(3)doi:10.3390/ijms21030736
# 105.  Ben Amara A, Gorvel L, Baulan K, et al. Placental macrophages are impaired in chorioamnionitis, an infectious pathology of the placenta. J Immunol. Dec 01 2013;191(11):5501-14. doi:10.4049/jimmunol.1300988

macrophage.lineage.markers.2 <- c("CD14", "ITGAM", "CSF1R", "ADGRE1", "CD80", "CD38", "GPR18", "FPR2", "CD86", "TNF", "IL12A", "MYC", "EGR2", "CD163", "MRC1", "CD209", "TGFB1", "IL10", "VEGFA", "CD163", "HLA-DRA", "CD86", "IL6")
macrophage.lineage.markers.2 <- unique(macrophage.lineage.markers.2)

feature.plot <- DotPlot(seurat.object, features = macrophage.lineage.markers.2)
png(file=paste0("macrophage.lineage.markers2-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Macrophage Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

# https://www.nature.com/articles/s41591-021-01329-2#Sec8
t.lineages <- c("CD4", "CD8A", "CCR7", "PTPRC", "SELL", "CD27", "CD38", "CD44", "CXCR5", "CD40LG", "CCR7", "FOXP3", "IKZF2", "TRGV9", "TRDV2", "TRAV10", "TRAJ18", "TRAV1-2", "CD3G", "FCGR3A", "NCAM1", "NCR1", "IFNG", "TBX21", "TNF", "GATA3", "IL4", "RORC", "IL17A", "IL17F", "IL21")
t.lineages <- unique(t.lineages)

# https://www.nature.com/articles/s41591-021-01329-2#Sec8
major.t <- c("CD4", "CD8A", "CCR7", "PTPRC", "SELL", "CD27", "CD38", "CD44", "CXCR5", "CD40LG", "CCR7", "FOXP3", "IKZF2")
gd.t<- c("TRGV9", "TRDV2")
MAIT.t <- c("TRAV10", "TRAJ18", "TRAV1-2")
NK.t <- c("CD3G", "FCGR3A", "NCAM1", "NCR1")
TH1.t <- c("IFNG", "TBX21", "TNF")
TH2.t <- c("GATA3", "IL4")
TH17.t <- c("RORC", "IL17A", "IL17F", "IL21")

seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
top5.heatmap <- DoHeatmap(seurat.object, features = t.lineages, raster = FALSE)
ggsave("t.lineages_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 11, height = 8, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = t.lineages, slot="counts", raster = FALSE)
ggsave("t.lineages_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 11, height = 8, units = "in")

top5.heatmap <- DoHeatmap(seurat.object, features = macrophage.lineage.markers.2, raster = FALSE)
ggsave("macrophage.lineage.markers2_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 11, height = 8, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = macrophage.lineage.markers.2, slot="counts", raster = FALSE)
ggsave("macrophage.lineage.markers2_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=t.lineages)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, slot="counts", assay="Spatial")
ggsave("t.lineages_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, assay="Spatial")
ggsave("t.lineages_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=macrophage.lineage.markers.2)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("macrophage.lineage.markers.2_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, assay="Spatial")
ggsave("macrophage.lineage.markers.2_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=lymphoid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("lymphoid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, assay="Spatial")
ggsave("lymphoid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=myeloid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("myeloid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, assay="Spatial")
ggsave("myeloid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


Idents(seurat.object) <- "seurat_clusters"
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
plot2<- SpatialDimPlot(seurat.object, images=image.list) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_cropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialDimPlot(seurat.object, images=image.list, crop=FALSE) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")


Idents(seurat.object) <- "Subtype"
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
plot2<- SpatialDimPlot(seurat.object, images=image.list) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_Subtype_cropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialDimPlot(seurat.object, images=image.list, crop=FALSE) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_Subtype_uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")





## http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_PI3K_AKT_MTOR_SIGNALING.html
inflammation <- list(c("ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "ADM", "ADORA2B", "ADRM1", "AHR", "APLNR", "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1", "BEST1", "BST2", "BTG2", "C3AR1", "C5AR1", "CALCRL", "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7", "CCRL2", "CD14", "CD40", "CD48", "CD55", "CD69", "CD70", "CD82", "CDKN1A", "CHST2", "CLEC5A", "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1", "CXCL10", "CXCL11", "CXCL6", "CXCL8", "CXCL9", "CXCR6", "CYBB", "DCBLD2", "EBI3", "EDN1", "EIF2AK2", "EMP3", "EREG", "F3", "FFAR2", "FPR1", "FZD5", "GABBR1", "GCH1", "GNA15", "GNAI3", "GP1BA", "GPC3", "GPR132", "GPR183", "HAS2", "HBEGF", "HIF1A", "HPN", "HRH1", "ICAM1", "ICAM4", "ICOSLG", "IFITM1", "IFNAR1", "IFNGR2", "IL10", "IL10RA", "IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5", "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2", "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV", "MEP1A", "MET", "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA", "NLRP3", "NMI", "NMUR1", "NOD2", "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2", "PCDH7", "PDE4B", "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2", "PTGER4", "PTGIR", "PTPRE", "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16", "RHOG", "RIPK2", "RNF144B", "ROS1", "RTP4", "SCARF1", "SCN1B", "SELE", "SELENOS", "SELL", "SEMA4D", "SERPINE1", "SGMS2", "SLAMF1", "SLC11A2", "SLC1A2", "SLC28A2", "SLC31A1", "SLC31A2", "SLC4A4", "SLC7A1", "SLC7A2", "SPHK1", "SRI", "STAB1", "TACR1", "TACR3", "TAPBP", "TIMP1", "TLR1", "TLR2", "TLR3", "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10", "TNFSF15", "TNFSF9", "TPBG", "VIP"))
seurat.object <- AddModuleScore(object = seurat.object, features = inflammation, ctrl = 5, name = 'inflammation.score')
p1 <- FeaturePlot(seurat.object, features = 'inflammation.score1')
ggsave("FeaturePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'inflammation.score1')
ggsave("RidgePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "inflammation.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "inflammation.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

ifn.alpha <- list(c("ADAR", "B2M", "BATF2", "BST2", "C1S", "CASP1", "CASP8", "CCRL2", "CD47", "CD74", "CMPK2", "CMTR1", "CNP", "CSF1", "CXCL10", "CXCL11", "DDX60", "DHX58", "EIF2AK2", "ELF1", "EPSTI1", "GBP2", "GBP4", "GMPR", "HELZ2", "HERC6", "HLA-C", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3", "IL15", "IL4R", "IL7", "IRF1", "IRF2", "IRF7", "IRF9", "ISG15", "ISG20", "LAMP3", "LAP3", "LGALS3BP", "LPAR6", "LY6E", "MOV10", "MVB12A", "MX1", "NCOA7", "NMI", "NUB1", "OAS1", "OASL", "OGFR", "PARP12", "PARP14", "PARP9", "PLSCR1", "PNPT1", "PROCR", "PSMA3", "PSMB8", "PSMB9", "PSME1", "PSME2", "RIPK2", "RNF31", "RSAD2", "RTP4", "SAMD9", "SAMD9L", "SELL", "SLC25A28", "SP110", "STAT2", "TAP1", "TDRD7", "TENT5A", "TMEM140", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TRIM5", "TXNIP", "UBA7", "UBE2L6", "USP18", "WARS1"))
seurat.object <- AddModuleScore(object = seurat.object, features = ifn.alpha, ctrl = 5, name = 'ifn.alpha.score')
p1 <- FeaturePlot(seurat.object, features = 'ifn.alpha.score1')
ggsave("FeaturePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'ifn.alpha.score1')
ggsave("RidgePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "ifn.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "ifn.alpha.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

ifn.beta <- list(c("ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM", "BST2", "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1", "CMPK2", "CMTR1", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX58", "DDX60", "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "GBP4", "GBP6", "GCH1", "GPR18", "GZMA", "HELZ2", "HERC6", "HIF1A", "HLA-A", "HLA-B", "HLA-DMA", "HLA-DQA1", "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IFNAR2", "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1", "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7", "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1", "METTL7B", "MT2A", "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12", "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP", "PNPT1", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1", "RIPK2", "RNF213", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP", "SERPING1", "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2", "SP110", "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4", "STAT1", "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP", "TDRD7", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8", "VCAM1", "WARS1", "XAF1", "XCL1", "ZBP1", "ZNFX1"))
seurat.object <- AddModuleScore(object = seurat.object, features = ifn.beta, ctrl = 5, name = 'ifn.beta.score')
p1 <- FeaturePlot(seurat.object, features = 'ifn.beta.score1')
ggsave("FeaturePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'ifn.beta.score1')
ggsave("RidgePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "ifn.beta.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "ifn.beta.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

apoptosis <- list(c("ADD1", "AIFM3", "ANKH", "ANXA1", "APP", "ATF3", "AVPR1A", "BAX", "BCAP31", "BCL10", "BCL2L1", "BCL2L10", "BCL2L11", "BCL2L2", "BGN", "BID", "BIK", "BIRC3", "BMF", "BMP2", "BNIP3L", "BRCA1", "BTG2", "BTG3", "CASP1", "CASP2", "CASP3", "CASP4", "CASP6", "CASP7", "CASP8", "CASP9", "CAV1", "CCNA1", "CCND1", "CCND2", "CD14", "CD2", "CD38", "CD44", "CD69", "CDC25B", "CDK2", "CDKN1A", "CDKN1B", "CFLAR", "CLU", "CREBBP", "CTH", "CTNNB1", "CYLD", "DAP", "DAP3", "DCN", "DDIT3", "DFFA", "DIABLO", "DNAJA1", "DNAJC3", "DNM1L", "DPYD", "EBP", "EGR3", "EMP1", "ENO2", "ERBB2", "ERBB3", "EREG", "ETF1", "F2", "F2R", "FAS", "FASLG", "FDXR", "FEZ1", "GADD45A", "GADD45B", "GCH1", "GNA15", "GPX1", "GPX3", "GPX4", "GSN", "GSR", "GSTM1", "GUCY2D", "H1-0", "HGF", "HMGB2", "HMOX1", "HSPB1", "IER3", "IFITM3", "IFNB1", "IFNGR1", "IGF2R", "IGFBP6", "IL18", "IL1A", "IL1B", "IL6", "IRF1", "ISG20", "JUN", "KRT18", "LEF1", "LGALS3", "LMNA", "LUM", "MADD", "MCL1", "MGMT", "MMP2", "NEDD9", "NEFH", "PAK1", "PDCD4", "PDGFRB", "PEA15", "PLAT", "PLCB2", "PLPPR4", "PMAIP1", "PPP2R5B", "PPP3R1", "PPT1", "PRF1", "PSEN1", "PSEN2", "PTK2", "RARA", "RELA", "RETSAT", "RHOB", "RHOT2", "RNASEL", "ROCK1", "SAT1", "SATB1", "SC5D", "SLC20A1", "SMAD7", "SOD1", "SOD2", "SPTAN1", "SQSTM1", "TAP1", "TGFB2", "TGFBR3", "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF12A", "TNFSF10", "TOP2A", "TSPO", "TXNIP", "VDAC2", "WEE1", "XIAP"))
seurat.object <- AddModuleScore(object = seurat.object, features = apoptosis, ctrl = 5, name = 'apoptosis.score')
p1 <- FeaturePlot(seurat.object, features = 'apoptosis.score1')
ggsave("FeaturePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'apoptosis.score1')
ggsave("RidgePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "apoptosis.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "apoptosis.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tnf.alpha <- list(c("ABCA1", "ACKR3", "AREG", "ATF3", "ATP2B1", "B4GALT1", "B4GALT5", "BCL2A1", "BCL3", "BCL6", "BHLHE40", "BIRC2", "BIRC3", "BMP2", "BTG1", "BTG2", "BTG3", "CCL2", "CCL20", "CCL4", "CCL5", "CCN1", "CCND1", "CCNL1", "CCRL2", "CD44", "CD69", "CD80", "CD83", "CDKN1A", "CEBPB", "CEBPD", "CFLAR", "CLCF1", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL11", "CXCL2", "CXCL3", "CXCL6", "DDX58", "DENND5A", "DNAJB4", "DRAM1", "DUSP1", "DUSP2", "DUSP4", "DUSP5", "EDN1", "EFNA1", "EGR1", "EGR2", "EGR3", "EHD1", "EIF1", "ETS2", "F2RL1", "F3", "FJX1", "FOS", "FOSB", "FOSL1", "FOSL2", "FUT4", "G0S2", "GADD45A", "GADD45B", "GCH1", "GEM", "GFPT2", "GPR183", "HBEGF", "HES1", "ICAM1", "ICOSLG", "ID2", "IER2", "IER3", "IER5", "IFIH1", "IFIT2", "IFNGR2", "IL12B", "IL15RA", "IL18", "IL1A", "IL1B", "IL23A", "IL6", "IL6ST", "IL7R", "INHBA", "IRF1", "IRS2", "JAG1", "JUN", "JUNB", "KDM6B", "KLF10", "KLF2", "KLF4", "KLF6", "KLF9", "KYNU", "LAMB3", "LDLR", "LIF", "LITAF", "MAFF", "MAP2K3", "MAP3K8", "MARCKS", "MCL1", "MSC", "MXD1", "MYC", "NAMPT", "NFAT5", "NFE2L2", "NFIL3", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NINJ1", "NR4A1", "NR4A2", "NR4A3", "OLR1", "PANX1", "PDE4B", "PDLIM5", "PER1", "PFKFB3", "PHLDA1", "PHLDA2", "PLAU", "PLAUR", "PLEK", "PLK2", "PLPP3", "PMEPA1", "PNRC1", "PPP1R15A", "PTGER4", "PTGS2", "PTPRE", "PTX3", "RCAN1", "REL", "RELA", "RELB", "RHOB", "RIPK2", "RNF19B", "SAT1", "SDC4", "SERPINB2", "SERPINB8", "SERPINE1", "SGK1", "SIK1", "SLC16A6", "SLC2A3", "SLC2A6", "SMAD3", "SNN", "SOCS3", "SOD2", "SPHK1", "SPSB1", "SQSTM1", "STAT5A", "TANK", "TAP1", "TGIF1", "TIPARP", "TLR2", "TNC", "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFRSF9", "TNFSF9", "TNIP1", "TNIP2", "TRAF1", "TRIB1", "TRIP10", "TSC22D1", "TUBB2A", "VEGFA", "YRDC", "ZBTB10", "ZC3H12A", "ZFP36"))
seurat.object <- AddModuleScore(object = seurat.object, features = tnf.alpha, ctrl = 5, name = 'tnf.alpha.score')
p1 <- FeaturePlot(seurat.object, features = 'tnf.alpha.score1')
ggsave("FeaturePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'tnf.alpha.score1')
ggsave("RidgePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "tnf.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "tnf.alpha.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")


hypoxia.list <- c("ACKR3", "ADM", "ADORA2B", "AK4", "AKAP12", "ALDOA", "ALDOB", "ALDOC", "AMPD3", "ANGPTL4", "ANKZF1", "ANXA2", "ATF3", "ATP7A", "B3GALT6", "B4GALNT2", "BCAN", "BCL2", "BGN", "BHLHE40", "BNIP3L", "BRS3", "BTG1", "CA12", "CASP6", "CAV1", "CAVIN1", "CAVIN3", "CCN1", "CCN2", "CCN5", "CCNG2", "CDKN1A", "CDKN1B", "CDKN1C", "CHST2", "CHST3", "CITED2", "COL5A1", "CP", "CSRP2", "CXCR4", "DCN", "DDIT3", "DDIT4", "DPYSL4", "DTNA", "DUSP1", "EDN2", "EFNA1", "EFNA3", "EGFR", "ENO1", "ENO2", "ENO3", "ERO1A", "ERRFI1", "ETS1", "EXT1", "F3", "FAM162A", "FBP1", "FOS", "FOSL2", "FOXO3", "GAA", "GALK1", "GAPDH", "GAPDHS", "GBE1", "GCK", "GCNT2", "GLRX", "GPC1", "GPC3", "GPC4", "GPI", "GRHPR", "GYS1", "HAS1", "HDLBP", "HEXA", "HK1", "HK2", "HMOX1", "HOXB9", "HS3ST1", "HSPA5", "IDS", "IER3", "IGFBP1", "IGFBP3", "IL6", "ILVBL", "INHA", "IRS2", "ISG20", "JMJD6", "JUN", "KDELR3", "KDM3A", "KIF5A", "KLF6", "KLF7", "KLHL24", "LALBA", "LARGE1", "LDHA", "LDHC", "LOX", "LXN", "MAFF", "MAP3K1", "MIF", "MT1E", "MT2A", "MXI1", "MYH9", "NAGK", "NCAN", "NDRG1", "NDST1", "NDST2", "NEDD4L", "NFIL3", "NOCT", "NR3C1", "P4HA1", "P4HA2", "PAM", "PCK1", "PDGFB", "PDK1", "PDK3", "PFKFB3", "PFKL", "PFKP", "PGAM2", "PGF", "PGK1", "PGM1", "PGM2", "PHKG1", "PIM1", "PKLR", "PKP1", "PLAC8", "PLAUR", "PLIN2", "PNRC1", "PPARGC1A", "PPFIA4", "PPP1R15A", "PPP1R3C", "PRDX5", "PRKCA", "PYGM", "RBPJ", "RORA", "RRAGD", "S100A4", "SAP30", "SCARB1", "SDC2", "SDC3", "SDC4", "SELENBP1", "SERPINE1", "SIAH2", "SLC25A1", "SLC2A1", "SLC2A3", "SLC2A5", "SLC37A4", "SLC6A6", "SRPX", "STBD1", "STC1", "STC2", "SULT2B1", "TES", "TGFB3", "TGFBI", "TGM2", "TIPARP", "TKTL1", "TMEM45A", "TNFAIP3", "TPBG", "TPD52", "TPI1", "TPST2", "UGP2", "VEGFA", "VHL", "VLDLR", "WSB1", "XPNPEP1", "ZFP36", "ZNF292")
seurat.object <- AddModuleScore(object = seurat.object, features = hypoxia.list, ctrl = 5, name = 'hypoxia.score')
p1 <- FeaturePlot(seurat.object, features = 'hypoxia.score1')
ggsave("FeaturePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'hypoxia.score1')
ggsave("RidgePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "hypoxia.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "hypoxia.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

complement.list <- c("ACTN2", "ADAM9", "ADRA2B", "AKAP10", "ANG", "ANXA5", "APOA4", "APOBEC3F", "APOBEC3G", "APOC1", "ATOX1", "BRPF3", "C1QA", "C1QC", "C1R", "C1S", "C2", "C3", "C4BPB", "C9", "CA2", "CALM1", "CALM3", "CASP1", "CASP10", "CASP3", "CASP4", "CASP5", "CASP7", "CASP9", "CBLB", "CCL5", "CD36", "CD40LG", "CD46", "CD55", "CD59", "CDA", "CDH13", "CDK5R1", "CEBPB", "CFB", "CFH", "CLU", "COL4A2", "CP", "CPM", "CPQ", "CR1", "CR2", "CSRP1", "CTSB", "CTSC", "CTSD", "CTSH", "CTSL", "CTSO", "CTSS", "CTSV", "CXCL1", "DGKG", "DGKH", "DOCK10", "DOCK4", "DOCK9", "DPP4", "DUSP5", "DUSP6", "DYRK2", "EHD1", "ERAP2", "F10", "F2", "F3", "F5", "F7", "F8", "FCER1G", "FCN1", "FDX1", "FN1", "FYN", "GATA3", "GCA", "GMFB", "GNAI2", "GNAI3", "GNB2", "GNB4", "GNG2", "GNGT2", "GP1BA", "GP9", "GPD2", "GRB2", "GZMA", "GZMB", "GZMK", "HNF4A", "HPCAL4", "HSPA1A", "HSPA5", "IL6", "IRF1", "IRF2", "IRF7", "ITGAM", "ITIH1", "JAK2", "KCNIP2", "KCNIP3", "KIF2A", "KLK1", "KLKB1", "KYNU", "L3MBTL4", "LAMP2", "LAP3", "LCK", "LCP2", "LGALS3", "LGMN", "LIPA", "LRP1", "LTA4H", "LTF", "LYN", "MAFF", "ME1", "MMP12", "MMP13", "MMP14", "MMP15", "MMP8", "MSRB1", "MT3", "NOTCH4", "OLR1", "PCLO", "PCSK9", "PDGFB", "PDP1", "PFN1", "PHEX", "PIK3CA", "PIK3CG", "PIK3R5", "PIM1", "PLA2G4A", "PLA2G7", "PLAT", "PLAUR", "PLEK", "PLG", "PLSCR1", "PPP2CB", "PPP4C", "PRCP", "PRDM4", "PREP", "PRKCD", "PRSS3", "PRSS36", "PSEN1", "PSMB9", "RABIF", "RAF1", "RASGRP1", "RBSN", "RCE1", "RHOG", "RNF4", "S100A12", "S100A13", "S100A9", "SCG3", "SERPINA1", "SERPINB2", "SERPINC1", "SERPINE1", "SERPING1", "SH2B3", "SIRT6", "SPOCK2", "SRC", "STX4", "TFPI2", "TIMP1", "TIMP2", "TMPRSS6", "TNFAIP3", "USP14", "USP15", "USP16", "USP8", "VCPIP1", "WAS", "XPNPEP1", "ZEB1", "ZFPM2")
seurat.object <- AddModuleScore(object = seurat.object, features = complement.list, ctrl = 5, name = 'complement.score')
p1 <- FeaturePlot(seurat.object, features = 'complement.score1')
ggsave("FeaturePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'complement.score1')
ggsave("RidgePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "complement.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "complement.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

angiogenesis.list <- c("APOH", "APP", "CCND2", "COL3A1", "COL5A2", "CXCL6", "FGFR1", "FSTL1", "ITGAV", "JAG1", "JAG2", "KCNJ8", "LPL", "LRPAP1", "LUM", "MSX1", "NRP1", "OLR1", "PDGFA", "PF4", "PGLYRP1", "POSTN", "PRG2", "PTK2", "S100A4", "SERPINA5", "SLCO2A1", "SPP1", "STC1", "THBD", "TIMP1", "TNFRSF21", "VAV2", "VCAN", "VEGFA", "VTN")
seurat.object <- AddModuleScore(object = seurat.object, features = angiogenesis.list, ctrl = 5, name = 'angiogenesis.score')
p1 <- FeaturePlot(seurat.object, features = 'angiogenesis.score1')
ggsave("FeaturePlot_angiogenesis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'angiogenesis.score1')
ggsave("RidgePlot_angiogenesis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "angiogenesis.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_angiogenesis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "angiogenesis.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_angiogenesis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")


il2 <- list(c("ABCB1", "ADAM19", "AGER", "AHCY", "AHNAK", "AHR", "ALCAM", "AMACR", "ANXA4", "APLP1", "ARL4A", "BATF", "BATF3", "BCL2", "BCL2L1", "BHLHE40", "BMP2", "BMPR2", "CA2", "CAPG", "CAPN3", "CASP3", "CCND2", "CCND3", "CCNE1", "CCR4", "CD44", "CD48", "CD79B", "CD81", "CD83", "CD86", "CDC42SE2", "CDC6", "CDCP1", "CDKN1C", "CISH", "CKAP4", "COCH", "COL6A1", "CSF1", "CSF2", "CST7", "CTLA4", "CTSZ", "CXCL10", "CYFIP1", "DCPS", "DENND5A", "DHRS3", "DRC1", "ECM1", "EEF1AKMT1", "EMP1", "ENO3", "ENPP1", "EOMES", "ETFBKMT", "ETV4", "F2RL2", "FAH", "FAM126B", "FGL2", "FLT3LG", "FURIN", "GABARAPL1", "GADD45B", "GALM", "GATA1", "GBP4", "GLIPR2", "GPR65", "GPR83", "GPX4", "GSTO1", "GUCY1B1", "HIPK2", "HK2", "HOPX", "HUWE1", "ICOS", "IFITM3", "IFNGR1", "IGF1R", "IGF2R", "IKZF2", "IKZF4", "IL10", "IL10RA", "IL13", "IL18R1", "IL1R2", "IL1RL1", "IL2RA", "IL2RB", "IL3RA", "IL4R", "IRF4", "IRF6", "IRF8", "ITGA6", "ITGAE", "ITGAV", "ITIH5", "KLF6", "LCLAT1", "LIF", "LRIG1", "LRRC8C", "LTB", "MAFF", "MAP3K8", "MAP6", "MAPKAPK2", "MUC1", "MXD1", "MYC", "MYO1C", "MYO1E", "NCOA3", "NCS1", "NDRG1", "NFIL3", "NFKBIZ", "NOP2", "NRP1", "NT5E", "ODC1", "P2RX4", "P4HA1", "PDCD2L", "PENK", "PHLDA1", "PHTF2", "PIM1", "PLAGL1", "PLEC", "PLIN2", "PLPP1", "PLSCR1", "PNP", "POU2F1", "PRAF2", "PRKCH", "PRNP", "PTCH1", "PTGER2", "PTH1R", "PTRH2", "PUS1", "RABGAP1L", "RGS16", "RHOB", "RHOH", "RNH1", "RORA", "RRAGD", "S100A1", "SCN9A", "SELL", "SELP", "SERPINB6", "SERPINC1", "SH3BGRL2", "SHE", "SLC1A5", "SLC29A2", "SLC2A3", "SLC39A8", "SMPDL3A", "SNX14", "SNX9", "SOCS1", "SOCS2", "SPP1", "SPRED2", "SPRY4", "ST3GAL4", "SWAP70", "SYNGR2", "SYT11", "TGM2", "TIAM1", "TLR7", "TNFRSF18", "TNFRSF1B", "TNFRSF21", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFSF10", "TNFSF11", "TRAF1", "TTC39B", "TWSG1", "UCK2", "UMPS", "WLS", "XBP1"))
seurat.object <- AddModuleScore(object = seurat.object, features = il2, ctrl = 5, name = 'il2.score')
p1 <- FeaturePlot(seurat.object, features = 'il2.score1')
ggsave("FeaturePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'il2.score1')
ggsave("RidgePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "il2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "il2.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tgfb <- list(c("ACVR1", "APC", "ARID4B", "BCAR3", "BMP2", "BMPR1A", "BMPR2", "CDH1", "CDK9", "CDKN1C", "CTNNB1", "ENG", "FKBP1A", "FNTA", "FURIN", "HDAC1", "HIPK2", "ID1", "ID2", "ID3", "IFNGR2", "JUNB", "KLF10", "LEFTY2", "LTBP2", "MAP3K7", "NCOR2", "NOG", "PMEPA1", "PPM1A", "PPP1CA", "PPP1R15A", "RAB31", "RHOA", "SERPINE1", "SKI", "SKIL", "SLC20A1", "SMAD1", "SMAD3", "SMAD6", "SMAD7", "SMURF1", "SMURF2", "SPTBN1", "TGFB1", "TGFBR1", "TGIF1", "THBS1", "TJP1", "TRIM33", "UBE2D3", "WWTR1", "XIAP"))
seurat.object <- AddModuleScore(object = seurat.object, features = tgfb, ctrl = 5, name = 'tgfb.score')
p1 <- FeaturePlot(seurat.object, features = 'tgfb.score1')
ggsave("FeaturePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'tgfb.score1')
ggsave("RidgePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "tgfb.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "tgfb.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

il6 <- list(c("A2M", "ACVR1B", "ACVRL1", "BAK1", "CBL", "CCL7", "CCR1", "CD14", "CD36", "CD38", "CD44", "CD9", "CNTFR", "CRLF2", "CSF1", "CSF2", "CSF2RA", "CSF2RB", "CSF3R", "CXCL1", "CXCL10", "CXCL11", "CXCL13", "CXCL3", "CXCL9", "DNTT", "EBI3", "FAS", "GRB2", "HAX1", "HMOX1", "IFNAR1", "IFNGR1", "IFNGR2", "IL10RB", "IL12RB1", "IL13RA1", "IL15RA", "IL17RA", "IL17RB", "IL18R1", "IL1B", "IL1R1", "IL1R2", "IL2RA", "IL2RG", "IL3RA", "IL4R", "IL6", "IL6ST", "IL7", "IL9R", "INHBE", "IRF1", "IRF9", "ITGA4", "ITGB3", "JUN", "LEPR", "LTB", "LTBR", "MAP3K8", "MYD88", "OSMR", "PDGFC", "PF4", "PIK3R5", "PIM1", "PLA2G2A", "PTPN1", "PTPN11", "PTPN2", "REG1A", "SOCS1", "SOCS3", "STAM2", "STAT1", "STAT2", "STAT3", "TGFB1", "TLR2", "TNF", "TNFRSF12A", "TNFRSF1A", "TNFRSF1B", "TNFRSF21", "TYK2"))
seurat.object <- AddModuleScore(object = seurat.object, features = il6, ctrl = 5, name = 'il6.score')
p1 <- FeaturePlot(seurat.object, features = 'il6.score1')
ggsave("FeaturePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'il6.score1')
ggsave("RidgePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "il6.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "il6.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

# http://www.gsea-msigdb.org/gsea/msigdb/cards/COATES_MACROPHAGE_M1_VS_M2_UP.html
M1 <- c("ABCA9", "ABHD1", "ACSS1", "ACYP2", "ADA", "AIF1", "ALAD", "ARFGEF3", "ARSB", "ATP6V0E2", "AXL", "BCKDHB", "BLOC1S6", "CADM1", "CAP1", "CAPN5", "CBX6", "CD59", "CFH", "CLBA1", "CNRIP1", "COLEC12", "COMT", "CRIM1", "CXCL14", "CXCR4", "DST", "DYNLT1", "EMC1", "ENO2", "FAM124A", "FAM135A", "FAM9B", "FGD2", "FILIP1L", "GALNT11", "GATM", "GDA", "GJA1", "GLO1", "GNB4", "HAUS2", "HDDC3", "HLA-DQA1", "HMGN3", "KCNJ10", "LAMA3", "LCORL", "LYPLAL1", "MAF", "MALAT1", "MARCKSL1", "MARCO", "MSR1", "NAT8L", "NRCAM", "OCEL1", "OGFRL1", "P2RY13", "PIANP", "PIK3AP1", "PLAAT3", "PLBD1", "PLXDC2", "PPP2R5C", "PTGER3", "RAB10", "RAPSN", "RASAL2", "RCBTB2", "RCN1", "RFX3", "RPL14", "SFI1", "SLC35A1", "SLC7A7", "SLCO2B1", "SRD5A3", "TGFBI", "TIFAB", "TM7SF3", "TOR3A", "TTC3", "TUBB2B", "TXNIP", "ZNF727")
seurat.object <- AddModuleScore(object = seurat.object, features = M1, ctrl = 5, name = 'M1.score')
p1 <- FeaturePlot(seurat.object, features = 'M1.score1')
ggsave("FeaturePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'M1.score1')
ggsave("RidgePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "M1.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "M1.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

M2<- c("ADAM8", "AK4", "AKR1E2", "ALDOC", "ANG", "ANKRD37", "ANOS1", "ARG1", "ARSK", "ATG4D", "ATP6V0D2", "BNIP3", "BST1", "C1QB", "C1QC", "C1QTNF12", "C5orf34", "CBLB", "CD24", "CD300C", "CD300LD", "CD5L", "CHIA", "COL20A1", "CRIPT", "CTSK", "CTSL", "CYTIP", "EFHD2", "EGLN3", "ERO1A", "F10", "F7", "FAM177A1", "FAM199X", "FAM241A", "FBLIM1", "FLRT2", "GASK1B", "GDF15", "GPRC5B", "HILPDA", "HLA-A", "HLA-B", "HS6ST1", "HSPA1A", "HSPA1B", "IGF2R", "IGHM", "IL18BP", "ITGB3", "LBP", "LIN7C", "LRRC27", "MCOLN3", "MFGE8", "MGST2", "MYO1F", "PADI4", "PDXDC1", "PINK1", "PKDCC", "PRDX2", "PROCR", "PTGER2", "QRICH1", "RGCC", "RIMBP2", "RPGRIP1", "RRAS2", "SCD", "SH2B2", "SLC2A1", "SNHG6", "SOAT1", "THBS1", "TJP2", "TLR1", "TMEM267", "TULP4", "UCHL1", "VEGFA", "XDH")
seurat.object <- AddModuleScore(object = seurat.object, features = M2, ctrl = 5, name = 'M2.score')
p1 <- FeaturePlot(seurat.object, features = 'M2.score1')
ggsave("FeaturePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'M2.score1')
ggsave("RidgePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "M2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "M2.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")


seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
save.image("human_spatial_data-macs-annotated_v1.RData")
# load("human_spatial_data-macs-annotated_v1.RData")
############################################################# ############################################################# #############################################################
############################################################# DE Analysis #############################################################
############################################################# ############################################################# #############################################################

seurat.object <- SCTransform(seurat.object, assay = "Spatial", verbose = T)
DefaultAssay(seurat.object) <- "Spatial"
seurat.object.var <- seurat.object@assays[["SCT"]]@var.features

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = seurat.object.var, assay="Spatial")

Idents(seurat.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-0.2lnFC.txt", sep="\t")

Idents(seurat.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object, features = seurat.object.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "integrated_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("integrated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("integrated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("integrated_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("integrated_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("integrated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

##################################################################################################################################################################
###################################################### Macrophage lineage annotations ######################################################
##################################################################################################################################################################
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

Idents(seurat.object) <- "seurat_clusters"
levels(seurat.object)
  ##                0   1       2     3     4     5     6     7     8     9     10    11      12
cluster.Type <- c("M2", "M2", "M0.M1", "M2a", "M0.M2", "M1", "M0.1", "M0.M1", "M0.M1", "M0", "DC", "M2.M1", "M1")
names(cluster.Type) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Type)
## add cluster Types as a metadata column
seurat.object$Polarization <- Idents(seurat.object)

## Re-order the clusters to group similar clusters together
Idents(seurat.object) <- "Polarization"
levels(seurat.object)
seurat.object2 <- subset(x = seurat.object, idents = c("M2", "M2", "M0.M1", "M2a", "M0.M2", "M1", "M0.1", "M0.M1", "M0.M1", "M0", "DC", "M2.M1", "M1"))
ncol(seurat.object2)
        # make sure you didn't lose any cells
        # 3621 cells
levels(seurat.object2)
levels(x=seurat.object2) <- c("M0", "M0.M1", "M1", "M0.M2", "M2", "DC")
levels(seurat.object2)
ncol(seurat.object)
    # 3621
seurat.object <- seurat.object2
rm(seurat.object2)
levels(seurat.object)

Idents(seurat.object) <- "seurat_clusters"
cluster.Type <- c("M2", "M2", "M0", "M2a", "M0", "M1", "M0", "M0", "M0", "M0", "M0", "M2", "M2a")
names(cluster.Type) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Type)
## add cluster Types as a metadata column
seurat.object$Lineage <- Idents(seurat.object)

## Re-order the clusters to group similar clusters together
Idents(seurat.object) <- "Lineage"
levels(seurat.object)
seurat.object2 <- subset(x = seurat.object, idents = c("M0", "M1", "M2", "M2a"))
ncol(seurat.object2)
        # make sure you didn't lose any cells
        # 3621 cells
levels(seurat.object2)
levels(x=seurat.object2) <- c("M0", "M1", "M2", "M2a")
levels(seurat.object2)
ncol(seurat.object)
    # 3621
seurat.object <- seurat.object2
rm(seurat.object2)
levels(seurat.object)

seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
save.image("human_spatial_data-macs-annotated_v1.RData")


Idents(seurat.object) <- "seurat_clusters"

Idents(seurat.object) <- "seurat_clusters"
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "Condition", label= "TRUE")
p1
p2 <- FeaturePlot(seurat.object, reduction = "umap", features = "percent.viral")
p2
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Lineage", label= "TRUE")
p3
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Polarization", label= "TRUE")
p4

umap.combined <- CombinePlots(plots = list(p1, p2, p3, p4))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("integrated_UMAP-combined_annotated.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")


image.list <- c("S01", "S03.1","S04.2", "S15.3")
p5 <- SpatialPlot(seurat.object, features = "SIGLEC1", images=image.list, ncol=2)
p5
ggsave("SIGLEC1_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S01", "S03.1","S04.2", "S15.3")
p5 <- SpatialPlot(seurat.object, features = "SIGLEC1", images=image.list, ncol=2, pt.size.factor = 4)
p5
ggsave("SIGLEC1_SpatialPlots-cropped-LARGESPOTS.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S01", "S03.1","S04.2", "S15.3")
p5 <- SpatialPlot(seurat.object, features = "SIGLEC1", images=image.list, ncol=2, pt.size.factor = 1.6, crop=FALSE)
p5
ggsave("SIGLEC1_SpatialPlots-cropped-med.SPOTS_uncropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


Idents(seurat.object) <- "Lineage"
image.list <- c("S01")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c("M0", "M1", "M2")), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.lineage_S01.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S03.1")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c("M0", "M1", "M2")), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.lineage_S03.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S04.2")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c("M0", "M1", "M2")), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.lineage_S04.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S15.3")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c("M0", "M1", "M2")), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.lineage_S15.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")

image.list <- c("S01")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c("M0", "M1", "M2")), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0, crop=FALSE)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.lineage_S01_uncropped.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S03.1")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c("M0", "M1", "M2")), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0, crop=FALSE)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.lineage_S03_uncropped.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S04.2")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c("M0", "M1", "M2")), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0, crop=FALSE)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.lineage_S04_uncropped.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S15.3")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c("M0", "M1", "M2")), facet.highlight = TRUE, ncol = 3, pt.size.factor = 4.0, crop=FALSE)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.lineage_S15_uncropped.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")

image.list <- c("S01")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c("M0", "M1", "M2")), facet.highlight = TRUE, ncol = 3, pt.size.factor = 6, crop=FALSE)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.lineage_S01_uncropped_XL.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S03.1")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c("M0", "M1", "M2")), facet.highlight = TRUE, ncol = 3, pt.size.factor = 6, crop=FALSE)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.lineage_S03_uncropped_XL.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S04.2")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c("M0", "M1", "M2")), facet.highlight = TRUE, ncol = 3, pt.size.factor = 6, crop=FALSE)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.lineage_S04_uncropped_XL.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S15.3")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c("M0", "M1", "M2")), facet.highlight = TRUE, ncol = 3, pt.size.factor = 6, crop=FALSE)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.lineage_S15_uncropped_XL.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")



seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
save.image("human_spatial_data-macs-annotated_v1.RData")

############################################################# ############################################################# #############################################################
############################################################# VISUALISATION OF TOP MARKERS #############################################################
############################################################# ############################################################# #############################################################
dir.create("Top.Markers")
setwd("Top.Markers")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S01", "S03.1","S04.2", "S15.3")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=2)
p5
ggsave("C3_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("S01", "S03.1","S04.2", "S15.3")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=2, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("C3_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("S01", "S03.1","S04.2", "S15.3")
p5 <- SpatialPlot(seurat.object, features = "KISS1", images=image.list, ncol=2)
p5
ggsave("KISS1_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("S01", "S03.1","S04.2", "S15.3")
p5 <- SpatialPlot(seurat.object, features = "KISS1", images=image.list, ncol=2, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("KISS1_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("S01", "S03.1","S04.2", "S15.3")
p5 <- SpatialPlot(seurat.object, features = "GDF15", images=image.list, ncol=2)
p5
ggsave("GDF15_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("S01", "S03.1","S04.2", "S15.3")
p5 <- SpatialPlot(seurat.object, features = "GDF15", images=image.list, ncol=2, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("GDF15_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


## SpatialPlots with all samples showing a particular genes # DAF is CD55, inhibits complement
image.list <- c("S01", "S03.1","S04.2", "S15.3")
marker.list <- c("C1QA", "C1QC", "C1R", "C1S", "C2", "C3", "CD55", "CD59")
p5 <- SpatialPlot(seurat.object, features = marker.list, images=image.list, ncol=2)
p5
ggsave("complement.select.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")
Idents(seurat.object) <- "orig.ident"
top5.heatmap <- DoHeatmap(seurat.object, features = marker.list, raster = FALSE)
ggsave("complement.select.markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 8, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = marker.list, slot="counts", raster = FALSE)
ggsave("complement.select.markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 8, units = "in")



marker.list <- c("C1QA", "C1QC", "C2", "C3", "CD55", "CFH", "F3", "KLK1", "SERPING1", "TIMP2", "EHD1")
p5 <- SpatialPlot(seurat.object, features = marker.list, images=image.list, ncol=2)
p5
ggsave("complement.top.pathway.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")

# http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_COMPLEMENT.html
feature.plot <- DotPlot(seurat.object, features = marker.list)
png(file=paste0("complement.select-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Complement Markers") + scale_y_discrete(name ="Condition")
dev.off()
marker.list <- c("ACTN2", "ADAM9", "ADRA2B", "AKAP10", "ANG", "ANXA5", "APOA4", "APOBEC3F", "APOBEC3G", "APOC1", "ATOX1", "BRPF3", "C1QA", "C1QC", "C1R", "C1S", "C2", "C3", "C4BPB", "C9", "CA2", "CALM1", "CALM3", "CASP1", "CASP10", "CASP3", "CASP4", "CASP5", "CASP7", "CASP9", "CBLB", "CCL5", "CD36", "CD40LG", "CD46", "CD55", "CD59", "CDA", "CDH13", "CDK5R1", "CEBPB", "CFB", "CFH", "CLU", "COL4A2", "CP", "CPM", "CPQ", "CR1", "CR2", "CSRP1", "CTSB", "CTSC", "CTSD", "CTSH", "CTSL", "CTSO", "CTSS", "CTSV", "CXCL1", "DGKG", "DGKH", "DOCK10", "DOCK4", "DOCK9", "DPP4", "DUSP5", "DUSP6", "DYRK2", "EHD1", "ERAP2", "F10", "F2", "F3", "F5", "F7", "F8", "FCER1G", "FCN1", "FDX1", "FN1", "FYN", "GATA3", "GCA", "GMFB", "GNAI2", "GNAI3", "GNB2", "GNB4", "GNG2", "GNGT2", "GP1BA", "GP9", "GPD2", "GRB2", "GZMA", "GZMB", "GZMK", "HNF4A", "HPCAL4", "HSPA1A", "HSPA5", "IL6", "IRF1", "IRF2", "IRF7", "ITGAM", "ITIH1", "JAK2", "KCNIP2", "KCNIP3", "KIF2A", "KLK1", "KLKB1", "KYNU", "L3MBTL4", "LAMP2", "LAP3", "LCK", "LCP2", "LGALS3", "LGMN", "LIPA", "LRP1", "LTA4H", "LTF", "LYN", "MAFF", "ME1", "MMP12", "MMP13", "MMP14", "MMP15", "MMP8", "MSRB1", "MT3", "NOTCH4", "OLR1", "PCLO", "PCSK9", "PDGFB", "PDP1", "PFN1", "PHEX", "PIK3CA", "PIK3CG", "PIK3R5", "PIM1", "PLA2G4A", "PLA2G7", "PLAT", "PLAUR", "PLEK", "PLG", "PLSCR1", "PPP2CB", "PPP4C", "PRCP", "PRDM4", "PREP", "PRKCD", "PRSS3", "PRSS36", "PSEN1", "PSMB9", "RABIF", "RAF1", "RASGRP1", "RBSN", "RCE1", "RHOG", "RNF4", "S100A12", "S100A13", "S100A9", "SCG3", "SERPINA1", "SERPINB2", "SERPINC1", "SERPINE1", "SERPING1", "SH2B3", "SIRT6", "SPOCK2", "SRC", "STX4", "TFPI2", "TIMP1", "TIMP2", "TMPRSS6", "TNFAIP3", "USP14", "USP15", "USP16", "USP8", "VCPIP1", "WAS", "XPNPEP1", "ZEB1", "ZFPM2")

# Warning message:
# In FetchData(object = object, vars = features, cells = cells) :
#  The following requested variables were not found (16 ): STX4, SERPINA1, S100A12, FCN1, ERAP2, CTSV, CR1, CD59, CASP5, CASP10, CA2, C4BPB, C1S, C1R, APOBEC3G, APOBEC3F

Idents(seurat.object) <- "orig.ident"

feature.plot <- DotPlot(seurat.object, features = marker.list)
png(file=paste0("complement.pathway-DotPlotorig.ident.png"),
                res=300, 
                width=7000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Complement Markers") + scale_y_discrete(name ="Condition")
dev.off()

top5.heatmap <- DoHeatmap(seurat.object, features = marker.list, raster = FALSE)
ggsave("orig.ident_complement.pathway_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 5, height = 20, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = marker.list, slot="counts", raster = FALSE)
ggsave("orig.ident_complement.pathway_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 5, height = 20, units = "in")



## SpatialPlots with all samples showing a particular gene
image.list <- c("S01", "S03.1","S04.2", "S15.3")
mock.markers <- c("LUM", "DKK1", "IGFBP1", "FN1", "PRG2")
p5 <- SpatialPlot(seurat.object, features = mock.markers, images=image.list, ncol=2)
p5
ggsave("Mock.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("S01", "S03.1","S04.2", "S15.3")
marker.list <- c("PSG7", "CSH2", "KISS1", "PSG1", "GDF15")
p5 <- SpatialPlot(seurat.object, features = marker.list, images=image.list, ncol=2)
p5
ggsave("SARS-CoV-2.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")


setwd("..")


############################################################# ############################################################# #############################################################
############################################################# Subset Analysis: symptomatic v asymptomatic #############################################################
############################################################# ############################################################# #############################################################
dir.create("ax.sx")
setwd("ax.sx")
    #dir.create("ax.sx.viral.sub")      ## Redo this using only the spots with viral transcripts
    #setwd("ax.sx.viral.sub")
            #dir.create("ax.sx.viral.sub.viralclust")      ## Redo this using only the spots with viral transcripts, and cluster using viral transcripts
             #setwd("ax.sx.viral.sub.viralclust")
Idents(seurat.object)
Idents(seurat.object) <- "Condition"
seurat.object2 <- subset(x = seurat.object, idents = c("SARS-CoV-2.ax", "SARS-CoV-2.sx"))
Idents(seurat.object2) <- "Condition"
  #   seurat.object3 <- subset(seurat.object2, percent.viral > 0)
  #   ncol(seurat.object3)    # 752
  #   seurat.object2 <- seurat.object3
  #   rm(seurat.object3)

ncol(seurat.object2)
  # 9446
seurat.object2 <- SCTransform(seurat.object2, assay = "Spatial", verbose = T)
DefaultAssay(seurat.object2) <- "SCT"
seurat.object2.var <- seurat.object2@assays[["SCT"]]@var.features

seurat.object2 <- ScaleData(seurat.object2, features = seurat.object2.var, assay = "SCT", vars.to.regress = c("CC.Difference", "nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"))
     #  seurat.object2 <- ScaleData(seurat.object2, features = seurat.object2.var, assay = "SCT")
          #  seurat.object2 <- ScaleData(seurat.object2, features = SARS.genes, assay = "SCT")

seurat.object2<- RunPCA(seurat.object2, features = seurat.object2.var, assay = "SCT", slot = "scale.data")
    # seurat.object2<- RunPCA(seurat.object2, features = SARS.genes, assay = "SCT", slot = "scale.data")

seurat.object2 <- FindNeighbors(seurat.object2, features = "seurat.object2.var", dims = 1:30)
    # seurat.object2 <- FindNeighbors(seurat.object2, features = "seurat.object2.var", dims = 1:5)
seurat.object2 <- FindClusters(seurat.object2, resolution = 0.6)
seurat.object2 <- RunUMAP(seurat.object2, dims=1:30)
    # seurat.object2 <- RunUMAP(seurat.object2, dims=1:5)
  ## clusters 0 thru 12
p1 <- DimPlot(seurat.object2, reduction = "umap", group.by = "orig.ident")
p1
p2 <- DimPlot(seurat.object2, reduction = "umap", group.by = "seurat_clusters", label= "TRUE")
p2
p3 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Type")
p3
p4 <- FeaturePlot(seurat.object2, reduction = "umap", features = "percent.viral")
p4
umap.combined <- CombinePlots(plots = list(p1, p2, p3, p4))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("integrated_UMAP.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object2) <- "Spatial"
seurat.object2 <- NormalizeData(seurat.object2, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object2 <- ScaleData(seurat.object2, features = all.genes, assay="Spatial")
seurat.object2 <- ScaleData(seurat.object2, features = all.genes, assay="SCT")

Idents(seurat.object2) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object2, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "integrated_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
de_markers <- FindAllMarkers(seurat.object2, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "integrated_DEGs_byclusters_pos-0.2lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
View(top5)
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, raster = FALSE)
ggsave("integrated_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("integrated_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("integrated_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
#ggsave("top2.marker-DotPlot.pdf", plot = feature.plot, device = "pdf", width = 16, height = 4, units = "in")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("integrated_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("integrated_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

##### Perform differential expression between Condition
Idents(seurat.object2) <- "Condition"
de_markers <- FindAllMarkers(seurat.object2, features = seurat.object2.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "Condition_DEGs_pos-0.693lnFC.txt", sep="\t")
de_markers <- FindAllMarkers(seurat.object2, features = seurat.object2.var, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "Condition_DEGs_pos-0.2lnFC.txt", sep="\t")

ggsave("Condition_top5_markers_logfc_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object2, features = top5$gene, slot="counts", raster = FALSE)
ggsave("Condition_top5_markers_counts_heatmaptop3k.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
View(top2)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object2, features = unique.top2)
png(file=paste0("Condition_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Condition")
dev.off()
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
plot2<- SpatialFeaturePlot(seurat.object2, features = "percent.viral", images=image.list) + patchwork::plot_layout(ncol = 4)
ggsave("infected-only_SpatialFeaturePlot_percent.viral.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

seurat.object <- ScaleData(seurat.object, features = SARS.genes, assay = "Spatial")


Idents(seurat.object) <- "Type"
feature.plot <- DotPlot(seurat.object, features = SARS.genes)
png(file=paste0("SARS.genes-DotPlot_Type.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_counts_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_logfc_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(seurat.object) <- "Condition"
feature.plot <- DotPlot(seurat.object, features = SARS.genes)
png(file=paste0("SARS.genes-DotPlot_Condition.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="Condition")
dev.off()
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_counts_Condition.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_logfc_Condition.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(seurat.object) <- "orig.ident"
feature.plot <- DotPlot(seurat.object, features = SARS.genes)
png(file=paste0("SARS.genes-DotPlot_Condition.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="orig.ident")
dev.off()
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = SARS.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes_seurat.clusters_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

rm(seurat.object2)

setwd("..")
############################################################# ############################################################# #############################################################
#############################################################  INTEGRATION WITH scRNA-seq DATASET   #############################################################
############################################################# ############################################################# #############################################################
############################################################# ############################################################# #############################################################
  ## remove all of the unnecessary lists
setwd("/home/ebarrozo/visium/results/seurat_human_v2")
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
# 

# save.image("human_spatial_data-integrated-annotated_v1.RData")


## see human_seurat_integrated_v1.R   for integration with tsang-saha_merged_v1.R dataset

