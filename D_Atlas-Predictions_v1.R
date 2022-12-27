## covid-visium_revisions_v1.R

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

## uses annotated and integrated visium object and re-analyzes based on new SARS-CoV-2 scheme
## integration required because clustering based on var.genes results in batch effect clustering
		## see viral.spatial.transcriptomics_v1.R and results in /home/ebarrozo/visium/results/seurat_human_v2/integration-v3/rcpa.integration/spatial/SCT_var.genes_UMAP
	## previous analysis was asymptomatic (ax) vs (sx)
	## new analysis 'case' will be positive vs negative for virus in the placenta tissue by Visium, RT-qPCR, RNAscope, and IHC probing for N and S
	## also analyzing by 'stage', which will separate presumed acute infections from subacute infections
	## lastly analyzing by 'viral', which will subset transcriptomes + for SARS-CoV-2 transcripts by Visium 
## Analyze data on AagaardLab3
# 



# https://satijalab.org/seurat/archive/v3.0/future_vignette.html
library(future)
plan()
# change the current plan to access parallelization
plan("multisession", workers = 40)
plan()
availableCores()
# Warning messages:
# 1: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead, explicitly specify either 'multisession' or 'multicore'. In the current R session, 'multiprocess' equals 'multisession'. 
# 2: In supportsMulticoreAndRStudio(...) :
 #  [ONE-TIME WARNING] Forked processing ('multicore') is not supported when running R from RStudio because it is considered unstable. For more details, how to control forked processing or not, and how to silence this warning in future R sessions, see ?parallelly::supportsMulticore
options(future.globals.maxSize = 220000000000)
	## now using 40 cores and max 220 gb ram

set.seed(seed=1)
library(dplyr)
library(Seurat)	 # Seurat_4.0.3
library(sctransform)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(glmGamPoi)	# BiocManager::install("glmGamPoi")
library(ggpubr)



setwd("/home/ebarrozo/visium/results/seurat_human_v2")
dir.create("revisions_v1")
setwd("revisions_v1")

## New folder for these revisions- see above
setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1")

dir.create("atlas")
setwd("atlas")

setwd("/home/ebarrozo/visium/results/seurat_human_v2/integration-v3/rcpa.integration")
load("human-rcpa-integrated-umap_v3.RData")
DefaultAssay(seurat.object) <- "RNA"
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap", assays=c("RNA"))

setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1/atlas")

library(ggpubr)
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "etal", cols = mypal3, raster=FALSE) + labs(title = NULL, color='et al.')
ggsave("UMAP_XL-etal.pdf", plot = p1, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Niche')
ggsave("UMAP-Type-Niche.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Annotation')
ggsave("UMAP-Type-Annotation.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")

p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Site", cols = mypal3, raster=FALSE) + labs(title = NULL, color='Site')
ggsave("UMAP-Site.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Platform", cols = mypal3, raster=FALSE) + labs(title = NULL, color='Platform')
ggsave("UMAP-Platform.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Condition", cols = mypal3, raster=FALSE) + labs(title = NULL, color='Condition')
ggsave("UMAP-Condition.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")

umap.combined <- p1 + p2 + p3 + p4
ggsave("integrated_UMAP-FINAL.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("UMAP_XL.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")


p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", label= "FALSE", raster=TRUE, cols = mypal3) + labs(title = NULL, color='Annotation')
ggsave("UMAP-Type-Annotation-rasterTRUE.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Site", cols = mypal3, raster=TRUE) + labs(title = NULL, color='Site')
ggsave("UMAP-Site-rasterTRUE.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Platform", cols = mypal3, raster=TRUE) + labs(title = NULL, color='Platform')
ggsave("UMAP-Platform-rasterTRUE.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Condition", cols = mypal3, raster=TRUE) + labs(title = NULL, color='Condition')
ggsave("UMAP-Condition-rasterTRUE.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")
umap.combined <- p1 + p2 + p3 + p4
ggsave("integrated_UMAP-FINAL-rasterTRUE.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("integrated_UMAP-FINAL-rasterTRUE_XL.pdf", plot = umap.combined, device = "pdf", width = 2, height = 12, units = "in")

Idents(seurat.object) <- "Platform"
levels(seurat.object)
ncol(seurat.object) # 291871
seurat.object.atlas <- subset(x = seurat.object, idents = c("scRNA-seq", "snRNA-seq"))
rm(seurat.object)
ncol(seurat.object.atlas)	#273944
		## This subset is giving fits, so go to human_seurat_integrated_v3.R and load everything but spatial and merge below

seurat.object.atlas <- merge(x = seurat.object, y = c(seurat.object.tsang, seurat.object.yang), merge.data = TRUE, project = "integrated")
rm(seurat.object)
rm(seurat.object.tsang)
rm(seurat.object.yang)
rm(object.11)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(umap.combined)
rm(C3m.genes)
rm(C4m.genes)
rm(cluster.Subtypes)
rm(cluster.Type)
rm(merged.var.combined)
rm(new.cluster.ids)
rm(new.metadata)
rm(new.row.names)
rm(pbmc.var.combined)
rm(new.col.names)
rm(seurat.object.var)
rm(t1)
rm(top.features)
rm(top2)
rm(top3k)
rm(tsang.pbmc.var.combined)
rm(unique.top2)
rm(var.combined)

DefaultAssay(seurat.object.atlas) <- "RNA"
seurat.object.atlas <- DietSeurat(seurat.object.atlas, counts=TRUE, data=TRUE, scale.data = FALSE, assays=c("RNA"))
ncol(seurat.object.atlas)
	# [1] 273944


setwd("/home/ebarrozo/visium/results/seurat_human_v2")
load("/home/ebarrozo/visium/results/seurat_human_v2/annotated/human_spatial_data-integrated-annotated_v1.RData")
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, assays=c("Spatial"))

## https://satijalab.org/seurat/articles/spatial_vignette.html
setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1/atlas")
## Add metadata below


## Run predictions for spatial transcriptome identities based on single-cell transcriptome annotations
Idents(seurat.object.atlas) <- "Type"
levels(seurat.object.atlas)
	##  [1] "Stromal"                "SYT"                    "EVT"                   
	#	 [4] "VEC"                    "Fibroblast"             "PVC"                   
	#	 [7] "VCT"                    "NK"                     "Endometrial_epithelium"
	#	[10] "DC"                     "Macrophage"             "RBC"                   
	#	[13] "CD8T"                   "B"                      "SmoothMuscle"          
	#	[16] "DC2"                    "Megakaryocyte"          "T"                     
	#	[19] "Endothelial" 

	# Idents(seurat.object.atlas) <- "etal"
# levels(seurat.object.atlas)
type.list <- 	levels(seurat.object.atlas)


## atlas w/o visium is 14.7gb and visium object is 3.5 gb


	# seurat.object.atlas <- SCTransform(seurat.object.atlas, ncells = 10000, verbose = FALSE)
seurat.object.atlas <- SCTransform(seurat.object.atlas, ncells = 3000, return.only.var.genes=TRUE, conserve.memory = TRUE, verbose = TRUE) 
seurat.object.atlas <- RunPCA(seurat.object.atlas, assay="SCT", verbose = TRUE)  
	# Warning: Adding image data that isn't associated with any assay present
		## ^ by design? 
		## try adding Assay=SCT
var.genes <- seurat.object.atlas@assays[["SCT"]]@var.features
seurat.object.atlas <- FindNeighbors(seurat.object.atlas, assay="SCT",  features = "var.genes", dims = 1:30)
seurat.object.atlas <- FindClusters(seurat.object.atlas, assay="SCT", resolution = 0.6)
RunUMAP(seurat.object.atlas, dims = 1:30)

seurat.object <- SCTransform(seurat.object, ncells = 3000, conserve.memory = TRUE, assay="Spatial", verbose = TRUE) 
seurat.object <- RunPCA(seurat.object, assay="SCT", verbose = TRUE)  
# var.genes <- seurat.object@assays[["SCT"]]@var.features
# seurat.object <- FindNeighbors(seurat.object, assay="SCT",  features = "var.genes", dims = 1:30)
# seurat.object <- FindClusters(seurat.object, assay="SCT", resolution = 0.6)
# seurat.object <-RunUMAP(seurat.object, dims = 1:30)

anchors <- FindTransferAnchors(reference = seurat.object.atlas, query = seurat.object, normalization.method = "SCT")
	## Warning: Adding image data that isn't associated with any assay present
predictions.assay <- TransferData(anchorset = anchors, refdata = seurat.object.atlas$Type, prediction.assay = TRUE, weight.reduction = seurat.object[["pca"]], dims = 1:30)

seurat.object[["Predictions"]] <- predictions.assay

DefaultAssay(seurat.object) <- "Predictions"
levels(seurat.object)
Idents(seurat.object) <- "orig.ident"

# Now we get prediction scores for each spot for each class. Of particular interest in the placenta are the macrophages and syncytial trophoblasts (SYTs). 
# Here we can distinguish between distinct regions for these subtypes, for example:
		# image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
image.list <- c("S01")

p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("SYT", "Macrophage"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
p1
ggsave("Predictions_SYT.Macrophage_S01.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S23b")
p1 <- SpatialFeaturePlot(seurat.object, images = image.list, features = c("SYT", "Macrophage"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
p1
ggsave("Predictions_SYT.Macrophage_S23b.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")


image.list <- c("S01")
seurat.object <- FindSpatiallyVariableFeatures(seurat.object, assay = "Predictions", selection.method = "markvariogram",
    features = rownames(seurat.object), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(seurat.object), 4)
	# "SmoothMuscle" "SYT"          "VCT"          "Fibroblast" 
p1 <- SpatialPlot(object = seurat.object, features = top.clusters, ncol = 2, images = image.list)
p1
ggsave("Predictions_top.clusters_S01.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

top.clusters <- SpatiallyVariableFeatures(seurat.object)
 # [1] "SmoothMuscle"           "SYT"                    "VCT"                   
 #[4] "Fibroblast"             "DC2"                    "B"                     
 #[7] "RBC"                    "max"                    "VEC"                   
#[10] "EVT"                    "Stromal"                "PVC"                   
# [13] "Macrophage"             "Endometrial-epithelium" "DC"

# rm max
top.clusters.1 <- top.clusters[-8]

image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
image.list <- c("S01")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S01.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S03.1")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S03.1.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S04.2")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S04.2.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S15.3")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S15.3.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")

image.list <- c("S16.4")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S16.4.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S17.5")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S17.5.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S18.6")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S18.6.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S19.7")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S19.7.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S20.8")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S20.8.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S21.9")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S21.9.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")

image.list <- c("S21.9")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S21.9.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S23a.11")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S23a.11.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S23b.12")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S23b.12.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S24.13")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S24.13.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S25.14")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S25.14.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S26.15")
p1<- SpatialPlot(seurat.object, features = top.clusters.1, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("Predictions_all_S26.15.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")






rownames((predictions.assay))
# [1] "EVT"                    "Stromal"                "CD8T"                  
# [4] "NK"                     "Macrophage"             "VCT"                   
 #[7] "Fibroblast"             "B"                      "SmoothMuscle"          
#[10] "DC2"                    "SYT"                    "Megakaryocyte"         
#[13] "VEC"                    "PVC"                    "Endometrial-epithelium"
#[16] "DC"                     "RBC"                    "T"                     
#[19] "Endothelial"            "max" 




Idents(seurat.object) <- "Predictions"
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
plot2<- SpatialDimPlot(seurat.object, images=image.list) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_cropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialDimPlot(seurat.object, images=image.list, crop=FALSE) + patchwork::plot_layout(ncol = 4)
ggsave("integrated_SpatialFeaturePlot_seurat_clusters_uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
 

# https://learn.gencore.bio.nyu.edu/seurat-integration-and-label-transfer/
seurat.object <- AddMetaData(object = seurat.object, metadata = Predictions)
seurat.object$prediction.match <- seurat.object$predicted.id == seurat.object$Type
table(seurat.object$prediction.match)
table(seurat.object$predicted.id)

# rm(predictions.assay)
rm(seurat.object.atlas)

rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(umap.combined)

# seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap", assays=c("Spatial", "Predictions"))

save.image("predictions_v1.RData")
	# load('predictions_v1.RData')	## 249 gb w/o removing atlas and predicitons.assay









setwd("/home/ebarrozo/visium/results/seurat_human_v2")
dir.create("revisions_v1")
setwd("revisions_v1")

## New folder for these revisions- see above
setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1")


## Check cell annotations
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", cols = mypal3, label=T, repel=T, raster=T) + labs(title = NULL, color='Niche')
p4   ## clusters are annotated
## Check seurat_clusters is still intact
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", cols = mypal3, label=T, repel=T, raster=T) + labs(title = NULL, color='Cluster')
p4   ## They are



## Add annotion data
new.metadata <- c("PV", "BP", "CAM", "PVBP", 
	"PVBP", "PVBP", "PVBP", "PVBP", 
	"PVBP", "PVBP", "PVBP", "PVBP", 
	"PVBP", "PVBP", "PVBP", "PVBP")
Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
	## [1] "S01"  "S03"  "S04"  "S15"  "S16"  "S17"  "S18"  "S19"  "S20"  "S21"  "S22"  "S23a" "S23b" "S24" "S25"  "S26" 
	## Any new annotation has to be in this order
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
## add cluster Subtypes as a metadata column
seurat.object$Site <- Idents(seurat.object)

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("Control", "Control", "Control", "SARS-CoV-2.Asymptomatic", "SARS-CoV-2.Symptomatic", "SARS-CoV-2.Asymptomatic", "SARS-CoV-2.Asymptomatic", "SARS-CoV-2.Asymptomatic", "SARS-CoV-2.Symptomatic", "SARS-CoV-2.Symptomatic", "SARS-CoV-2.Symptomatic", "SARS-CoV-2.Symptomatic", "SARS-CoV-2.Symptomatic", "Control","Control","Control")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Condition <- Idents(seurat.object)
Idents(seurat.object) <- "Condition"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("F", "F", "F", "F", 
	"M", "F", "M", "M", 
	"F", "F", "M", "F", 
	"F", "M","M","M")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$FetalSex <- Idents(seurat.object)
Idents(seurat.object) <- "FetalSex"

## New annotations 8/3/2022
Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
 new.metadata <- c("Control", "Control", "Control", "Respiratory.Asymptomatic", "Placenta.Subacute", "Respiratory.Asymptomatic", "Respiratory.Asymptomatic", "Respiratory.Asymptomatic", "Respiratory.Symptomatic", "Respiratory.Symptomatic", "Placenta.Subacute", "Placenta.Acute", "Placenta.Acute", "Control","Control","Control")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Class <- Idents(seurat.object)
Idents(seurat.object) <- "Class"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
 new.metadata <- c("Control", "Control", "Control", "NotDetected", "Detected", "Detected", "Detected", "NotDetected", "NotDetected", "NotDetected", "Detected", "Detected", "Detected", "Control","Control","Control")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Orthogonal.Detection <- Idents(seurat.object)
Idents(seurat.object) <- "Orthogonal.Detection"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
 new.metadata <- c("Control", "Control", "Control", "NotDetected", "Detected", "NotDetected", "NotDetected", "NotDetected", "NotDetected", "NotDetected", "Detected", "Detected", "Detected", "Control","Control","Control")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Visium.Detection <- Idents(seurat.object)
Idents(seurat.object) <- "Visium.Detection"



####################################################################################################
####################################################################################################
####################################################################################################

## Sexual_dimorphisms in subsets of visium data



## since S23 was acute and IUFD, we analyzed all other SARS-CoV-2 samples independently for sexual dimorphism differences. 
Idents(seurat.object) <- "Class"
levels(seurat.object)   # "Control", "Respiratory.Asymptomatic", "Placenta.Subacute", Placenta.Acute , Respiratory.Symptomatic"
seurat.object.detected <- subset(x = seurat.object, idents = c("Placenta.Subacute", "Respiratory.Symptomatic", "Respiratory.Asymptomatic"))


dir.create("DE_FetalSex_S15-S22")
setwd("DE_FetalSex_S15-S22")
seurat.object.x <- seurat.object.detected
Idents(seurat.object.x) <- "FetalSex"
de_markers <- FindAllMarkers(seurat.object.x, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_FetalSex_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object.x, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
feature.plot
png(file=paste0("FetalSex_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
png(file=paste0("FetalSex_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
png(file=paste0("FetalSex_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object.x, features = top2$gene, raster = FALSE, angle=90)
ggsave("FetalSex_top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.x, features = top2$gene, raster = FALSE, slot="counts", angle=90)
ggsave("FetalSex_top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = FALSE, assay="Spatial", slot="counts")
ggsave("top30_counts_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("top30_fold-change_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.Top10_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("FetalSex.Top10.2_FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("FetalSex.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")

de_markers <- FindAllMarkers(seurat.object.x, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_FetalSex_0.2lnFC.txt", sep="\t")
setwd("..")

## Slim down the seurat.object.xs to save space. Don't need to keep the scale.data
seurat.object.x <- DietSeurat(seurat.object.x, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

rm(seurat.object.x)
rm(seurat.object.detected)

Idents(seurat.object) <- "Condition"
seurat.object.ax <- subset(x = seurat.object, idents = c("SARS-CoV-2.Asymptomatic"))
seurat.object.sx <- subset(x = seurat.object, idents = c("SARS-CoV-2.Symptomatic"))


Idents(seurat.object) <- "Visium.Detection"
seurat.object.detected <- subset(x = seurat.object, idents = c("Detected"))
seurat.object.notdetected <- subset(x = seurat.object, idents = c("NotDetected"))


dir.create("DE_FetalSex_Detected")
setwd("DE_FetalSex_Detected")
seurat.object.x <- seurat.object.detected
Idents(seurat.object.x) <- "FetalSex"
de_markers <- FindAllMarkers(seurat.object.x, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_FetalSex_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object.x, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
feature.plot
png(file=paste0("FetalSex_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
png(file=paste0("FetalSex_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
png(file=paste0("FetalSex_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object.x, features = top2$gene, raster = FALSE, angle=90)
ggsave("FetalSex_top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.x, features = top2$gene, raster = FALSE, slot="counts", angle=90)
ggsave("FetalSex_top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.Top10_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("FetalSex.Top10.2_FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("FetalSex.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = FALSE, assay="Spatial", slot="counts")
ggsave("top30_counts_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("top30_fold-change_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


de_markers <- FindAllMarkers(seurat.object.x, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_FetalSex_0.2lnFC.txt", sep="\t")
setwd("..")

## Slim down the seurat.object.xs to save space. Don't need to keep the scale.data
seurat.object.x <- DietSeurat(seurat.object.x, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")


dir.create("DE_FetalSex_NotDetected")
setwd("DE_FetalSex_NotDetected")
seurat.object.x <- seurat.object.notdetected
Idents(seurat.object.x) <- "FetalSex"
de_markers <- FindAllMarkers(seurat.object.x, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_FetalSex_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object.x, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
feature.plot
png(file=paste0("FetalSex_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
png(file=paste0("FetalSex_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
png(file=paste0("FetalSex_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object.x, features = top2$gene, raster = FALSE, angle=90)
ggsave("FetalSex_top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.x, features = top2$gene, raster = FALSE, slot="counts", angle=90)
ggsave("FetalSex_top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.Top10_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("FetalSex.Top10.2_FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("FetalSex.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = FALSE, assay="Spatial", slot="counts")
ggsave("top30_counts_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("top30_fold-change_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


de_markers <- FindAllMarkers(seurat.object.x, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_FetalSex_0.2lnFC.txt", sep="\t")
setwd("..")

## Slim down the seurat.object.xs to save space. Don't need to keep the scale.data
seurat.object.x <- DietSeurat(seurat.object.x, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")


dir.create("DE_FetalSex_ax")
setwd("DE_FetalSex_ax")
seurat.object.x <- seurat.object.ax
Idents(seurat.object.x) <- "FetalSex"
de_markers <- FindAllMarkers(seurat.object.x, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_FetalSex_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object.x, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
feature.plot
png(file=paste0("FetalSex_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
png(file=paste0("FetalSex_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
png(file=paste0("FetalSex_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object.x, features = top2$gene, raster = FALSE, angle=90)
ggsave("FetalSex_top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.x, features = top2$gene, raster = FALSE, slot="counts", angle=90)
ggsave("FetalSex_top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.Top10_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("FetalSex.Top10.2_FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("FetalSex.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = FALSE, assay="Spatial", slot="counts")
ggsave("top30_counts_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("top30_fold-change_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


de_markers <- FindAllMarkers(seurat.object.x, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_FetalSex_0.2lnFC.txt", sep="\t")
setwd("..")

## Slim down the seurat.object.xs to save space. Don't need to keep the scale.data
seurat.object.x <- DietSeurat(seurat.object.x, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")





dir.create("DE_FetalSex_sx")
setwd("DE_FetalSex_sx")
seurat.object.x <- seurat.object.sx
Idents(seurat.object.x) <- "FetalSex"
de_markers <- FindAllMarkers(seurat.object.x, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_FetalSex_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object.x, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
feature.plot
png(file=paste0("FetalSex_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
png(file=paste0("FetalSex_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
png(file=paste0("FetalSex_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object.x, features = top2$gene, raster = FALSE, angle=90)
ggsave("FetalSex_top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.x, features = top2$gene, raster = FALSE, slot="counts", angle=90)
ggsave("FetalSex_top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.Top10_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("FetalSex.Top10.2_FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("FetalSex.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = FALSE, assay="Spatial", slot="counts")
ggsave("top30_counts_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("top30_fold-change_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


de_markers <- FindAllMarkers(seurat.object.x, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_FetalSex_0.2lnFC.txt", sep="\t")
setwd("..")

## Slim down the seurat.object.xs to save space. Don't need to keep the scale.data
seurat.object.x <- DietSeurat(seurat.object.x, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")


rm(seurat.object.x)
rm(seurat.object.ax)
rm(seurat.object.sx)
rm(seurat.object.detected)
rm(seurat.object.notdetected)
####################################################################################################
####################################################################################################
####################################################################################################



## Save UMAPs
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Niche')
p2
ggsave("UMAP-Type-Niche.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Subtype",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Niche')
p2
ggsave("UMAP-Subype-Niche.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Annotation')
p2
ggsave("UMAP-Type-Annotation.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Subtype",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Annotation')
p2
ggsave("UMAP-Subype-Annotation.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Class",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Class')
p2
ggsave("UMAP-Class.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
ggsave("UMAP-Class_XL.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")

p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Orthogonal.Detection",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Orthogonal.Detection')
p2
ggsave("UMAP-Orthogonal.Detection.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Visium.Detection",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Visium.Detection')
p2
ggsave("UMAP-Visium.Detection.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Condition",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Condition')
p2
ggsave("UMAP-Condition.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
ggsave("UMAP-Condition_XL.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")




####################################################################################################
## Subset only spots with viral RNA
viral.object <- subset(seurat.object, subset = percent.viral > 0.000001)
ncol(viral.object)
	# 752 spots
write.table(viral.object@assays[["Spatial"]]@counts, "viral.counts.txt", sep="\t")
# Save a table with cell counts per cluster
Idents(viral.object) <- "orig.ident"
unfiltered.count <- table(Idents(viral.object))
write.table(unfiltered.count, "viral.object_orig.ident_cellcounts.txt", sep="\t")
	## S16 1, S22 1, S23a 15, S23b 735
## Subset an object with only virus genes
viral.object.viral.genes <- subset(viral.object, features = SARS.genes)
write.csv(viral.object.viral.genes@assays[["Spatial"]]@counts, file = "orig.ident_viral.counts.csv")
Idents(viral.object.viral.genes) <- "Class"
cluster.averages <- AverageExpression(viral.object.viral.genes, assays="Spatial", slot="counts", return.seurat = TRUE, features=SARS.genes)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("Viral.only-Class.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("Viral.only-Class.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
rm(viral.object.viral.genes)
rm(viral.object)

####################################################################################################
setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1")

####################################################################################################
######################### Perform differential expression between X   ########################################
########################################################################################################################
#### Start doing some Differential Expression ; Instead of DoScale and Normalize, Visium requires this below. Takes awhile, so run once before all DE, then diet before saving.
seurat.object <- SCTransform(seurat.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes = FALSE, verbose = T)
DefaultAssay(seurat.object) <- "SCT"
var.genes <- seurat.object@assays[["SCT"]]@var.features
all.genes <- rownames(seurat.object@assays[["Spatial"]])

## Run this before iterations then diet 
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes, assay="Spatial")

####################################################################################################
######################### Perform differential expression between Class   ########################################
########################################################################################################################
dir.create("DE_Class")
setwd("DE_Class")
Idents(seurat.object) <- "Class"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_Class_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
feature.plot
png(file=paste0("Class_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Class")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Class_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Class")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Class_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Class")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = top2$gene, raster = FALSE, angle=90)
ggsave("Class_top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object, features = top2$gene, raster = FALSE, slot="counts", angle=90)
ggsave("Class_top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("Class.Top10_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("Class.Top10.2_FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("Class.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("Class.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")

de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_Class_0.2lnFC.txt", sep="\t")
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
feature.plot
png(file=paste0("Class_top2.marker-DotPlot0.2lnFC.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Class")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Class_top3.marker-DotPlot0.2lnFC.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Class")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Class_top5.marker-DotPlot0.2lnFC.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Class")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = top2$gene, raster = FALSE, angle=90)
ggsave("Class_top5_markers_logfc_heatmap0.2lnFC.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object, features = top2$gene, raster = FALSE, slot="counts", angle=90)
ggsave("Class_top5_markers_counts_heatmap0.2lnFC.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("Class.Top10_Counts-MeanCounts_heatmap-0.2lnFC.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("Class.Top10.2_FoldChange-MeanCounts_heatmap-0.2lnFC.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
setwd("..")
####################################################################################################
######################### Perform differential expression between Orthogonal.Detection   ########################################
########################################################################################################################
dir.create("DE_Orthogonal.Detection")
setwd("DE_Orthogonal.Detection")
Idents(seurat.object) <- "Orthogonal.Detection"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_Orthogonal.Detection_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
feature.plot
png(file=paste0("Orthogonal.Detection_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Orthogonal.Detection")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Orthogonal.Detection_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Orthogonal.Detection")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Orthogonal.Detection_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Orthogonal.Detection")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = top2$gene, raster = FALSE, angle=90)
ggsave("Orthogonal.Detection_top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object, features = top2$gene, raster = FALSE, slot="counts", angle=90)
ggsave("Orthogonal.Detection_top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("Orthogonal.Detection.Top10_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("Orthogonal.Detection.Top10.2_FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("Orthogonal.Detection.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("Orthogonal.Detection.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")

de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_Orthogonal.Detection_0.2lnFC.txt", sep="\t")
setwd("..")
####################################################################################################
######################### Perform differential expression between Visium.Detection   ########################################
########################################################################################################################
dir.create("DE_Visium.Detection")
setwd("DE_Visium.Detection")
Idents(seurat.object) <- "Visium.Detection"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_Visium.Detection_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
feature.plot
png(file=paste0("Visium.Detection_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Visium.Detection")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Visium.Detection_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Visium.Detection")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Visium.Detection_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Visium.Detection")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = top2$gene, raster = FALSE, angle=90)
ggsave("Visium.Detection_top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object, features = top2$gene, raster = FALSE, slot="counts", angle=90)
ggsave("Visium.Detection_top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("Visium.Detection.Top10_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("Visium.Detection.Top10.2_FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("Visium.Detection.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("Visium.Detection.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")

de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_Visium.Detection_0.2lnFC.txt", sep="\t")
setwd("..")
####################################################################################################
######################### Perform differential expression between Condition   ########################################
########################################################################################################################
dir.create("DE_Condition")
setwd("DE_Condition")
Idents(seurat.object) <- "Condition"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_Condition_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
feature.plot
png(file=paste0("Condition_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Condition")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Condition_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Condition")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Condition_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Condition")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = top2$gene, raster = FALSE, angle=90)
ggsave("Condition_top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object, features = top2$gene, raster = FALSE, slot="counts", angle=90)
ggsave("Condition_top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("Condition.Top10_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("Condition.Top10.2_FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("Condition.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("Condition.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")

de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_Condition_0.2lnFC.txt", sep="\t")
setwd("..")
####################################################################################################
######################### Perform differential expression between FetalSex   ########################################
########################################################################################################################
dir.create("DE_FetalSex")
setwd("DE_FetalSex")
Idents(seurat.object) <- "FetalSex"
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes, assay="Spatial")

dimorphism.list <- c("FCGRT", "FCGR1A", "FCGR3B","IFI6", "CXCL10", "OAS1", "CCL2", "MX1", "IL10", "IFNA1", "IFNA2" ,"IFNG", "CCL4", "CD163")
Idents(seurat.object) <- "FetalSex"
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=dimorphism.list)
write.csv(cluster.averages@assays[["Spatial"]]@counts, file = "dimorphism.list_meancounts.csv")
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = dimorphism.list, raster = FALSE, slot="counts", assay="Spatial")
ggsave("dimorphism.list_mc_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = dimorphism.list, raster = FALSE, assay="Spatial")
ggsave("dimorphism.list_FC_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = dimorphism.list, raster = FALSE, assay="Spatial", slot="counts")
ggsave("dimorphism.list_counts_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = dimorphism.list, raster = FALSE, assay="Spatial")
ggsave("dimorphism.list_fold-change_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

de_markers <- FindAllMarkers(seurat.object, features = dimorphism.list, assay = "Spatial", only.pos = F, min.pct = 0, logfc.threshold =  0)
write.table(de_markers, "DEGs_by_FetalSex_dimorphism.list.txt", sep="\t")

de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_FetalSex_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
feature.plot
png(file=paste0("FetalSex_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("FetalSex_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("FetalSex_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object, features = top2$gene, raster = FALSE, angle=90)
ggsave("FetalSex_top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object, features = top2$gene, raster = FALSE, slot="counts", angle=90)
ggsave("FetalSex_top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in", size=3) + NoLegend()
## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.Top10_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("FetalSex.Top10.2_FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("FetalSex.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")

de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_FetalSex_0.2lnFC.txt", sep="\t")
setwd("..")

## Slim down the seurat.objects to save space. Don't need to keep the scale.data
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")


## Take the top 50 most upregulated genes in each class, based on fold-change and significance filter
acute.markers <- c("ORF10", "ORF1AB-1", "N", "TIMP1", "MT2A", "PRRG3", "CXCL10", "GPM6A", "ISG15", "CXCL8", "SOD2", "SLC2A3", "NUPR1", "FTH1", "SERPINE1", "MTRNR2L12", "MTRNR2L8", "LDHA", "NAMPT", "HK2", "IFIT3", "C15ORF48", "HLA-C", "BHLHE40", "GBP1", "IL1RL1", "UBC", "IL32", "MIF", "SLC25A37", "ENO1", "DDIT3", "IFITM3", "AL627171.2", "GAPDH", "STC1", "S100A9", "BNIP3", "JUNB", "FAM162A", "B2M", "HILPDA", "PNRC1", "GBP2", "GPI", "ZFAS1", "OPTN", "LYZ", "STAT1", "IRF1")
subacute.markers <- c("PSG7", "INHBA", "PAPPA2", "CRH", "CGA", "GDF15", "PSG6", "PSG3", "KISS1", "CGB3", "CSH1", "EBI3", "PSG11", "GH2", "PSG1", "HTRA1", "TFPI2", "PSG9", "PSG2", "ADAM12", "CYP19A1", "PSG4", "PSG5", "HPGD", "MMP11", "CCSAP", "SPTLC3", "FLT1", "EFHD1", "MAN1C1", "LINC00967", "TIMP2", "HMGB3", "ERVH48-1", "SDC1", "ERRFI1", "TFRC", "SLCO2A1", "PLIN2", "HSD3B1", "TFPI", "PSG8", "NECTIN4", "PAPPA", "SEMA3B", "KRT19", "TAC3", "SYDE1", "HOPX", "SLC2A1")
respiratory.symptomatic <- c("CSH2", "TFPI2", "PSG2", "GH2", "PSG5", "PSG4", "PSG1", "SERPINB2", "CGA", "CYP19A1", "ALPP", "SDC1", "HOPX", "PAPPA", "STS", "FBLN1", "PAPPA2", "XIST", "KISS1")
respiratory.asymptomatic <- c("CCN2", "XIST", "APOLD1", "FCGR2B", "PSG7", "FCGBP", "THBS1")
inf.notdetected <- union(respiratory.asymptomatic, respiratory.symptomatic)
inf.notdetected <- unique(inf.notdetected)
inf.detected <- union(acute.markers, subacute.markers)
inf.detected <- unique(inf.detected)
control.markers <- c("IGFBP1", "PRL", "DKK1", "LUM", "DCN", "FOS", "CHRDL1", "RBP1", "TIMP3", "CD248", "LGALS3", "GAS1", "HSD11B1", "IGFBP2", "ISLR", "MT-CO3", "FN1", "ABI3BP", "AOC1", "PRG2", "CRLF1", "CRIP1", "SCARA5", "MT-CO2", "MT-CYB", "MT-ND5", "MT-ND4", "IGFBP6", "MT-CO1", "PTGDS", "MT-ND1", "MFAP4", "FABP4", "S100A4", "IGFBP3", "MT-ATP6", "HBB", "MT-ND4L", "MT-ND3", "S100A6", "HLA-G", "NPW", "SERPINA3", "APOD", "CRISPLD2", "C3", "FOSB", "CFD", "MT-ND2", "CORO6")

####################################################################################################
## Subset only spots with viral RNA
viral.object <- subset(seurat.object, subset = percent.viral > 0.000001)
ncol(viral.object)
	# 752 spots
write.table(viral.object@assays[["Spatial"]]@counts, "viral.counts.txt", sep="\t")
# Save a table with cell counts per cluster
Idents(viral.object) <- "orig.ident"
unfiltered.count <- table(Idents(viral.object))
write.table(unfiltered.count, "viral.object_orig.ident_cellcounts.txt", sep="\t")
	## S16 1, S22 1, S23a 15, S23b 735
## Subset an object with only virus genes
viral.object.viral.genes <- subset(viral.object, features = SARS.genes)
write.csv(viral.object.viral.genes@assays[["Spatial"]]@counts, file = "orig.ident_viral.counts.csv")
Idents(viral.object.viral.genes) <- "Class"
cluster.averages <- AverageExpression(viral.object.viral.genes, assays="Spatial", slot="counts", return.seurat = TRUE, features=SARS.genes)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("Viral.only-Class.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("Viral.only-Class.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
rm(viral.object.viral.genes)
rm(viral.object)

####################################################################################################
setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1")

dir.create("viral.only")
setwd("viral.only")
viral.object <- subset(seurat.object, subset = percent.viral > 0.000001)


ncol(viral.object)
	# 752 spots
## Run this before iterations then diet 
DefaultAssay(viral.object) <- "Spatial"
viral.object <- NormalizeData(viral.object, normalization.method = "LogNormalize", scale.factor = 10000)
viral.object <- ScaleData(viral.object, features = all.genes, assay="Spatial")
## make essential figs with the cluster Subtype
## Save UMAPs
p2 <- DimPlot(viral.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(viral.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Niche')
p2
ggsave("UMAP-Type-Niche.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(viral.object, reduction = "umap", group.by = "Subtype",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Niche')
p2
ggsave("UMAP-Subype-Niche.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(viral.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Annotation')
p2
ggsave("UMAP-Type-Annotation.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(viral.object, reduction = "umap", group.by = "Subtype",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Annotation')
p2
ggsave("UMAP-Subype-Annotation.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(viral.object, reduction = "umap", group.by = "Class",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Class')
p2
ggsave("UMAP-Class.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
ggsave("UMAP-Class_XL.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")

p2 <- DimPlot(viral.object, reduction = "umap", group.by = "Orthogonal.Detection",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Orthogonal.Detection')
p2
ggsave("UMAP-Orthogonal.Detection.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(viral.object, reduction = "umap", group.by = "Visium.Detection",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Visium.Detection')
p2
ggsave("UMAP-Visium.Detection.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(viral.object, reduction = "umap", group.by = "Condition",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Condition')
p2
ggsave("UMAP-Condition.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
ggsave("UMAP-Condition_XL.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")

Idents(viral.object) <- "Type"
feature.plot <- DotPlot(viral.object, features = SARS.genes)
png(file=paste0("SARS.genes-Type-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="Niche")
dev.off()
Idents(viral.object) <- "orig.ident"
library.averages.heatmap <- DoHeatmap(viral.object, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes_orig.ident_heatmap_counts.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = SARS.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes_orig.ident_heatmap_logfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(viral.object) <- "Type"
feature.plot <- DotPlot(viral.object, features = acute.markers)
png(file=paste0("acute.markers-Type-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="Niche")
dev.off()
Idents(viral.object) <- "orig.ident"
library.averages.heatmap <- DoHeatmap(viral.object, features = acute.markers, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("acute.markers_orig.ident_heatmap_counts.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = acute.markers, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("acute.markers_orig.ident_heatmap_logfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(viral.object) <- "Type"
feature.plot <- DotPlot(viral.object, features = subacute.markers)
png(file=paste0("subacute.markers-Type-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="Niche")
dev.off()
Idents(viral.object) <- "orig.ident"
library.averages.heatmap <- DoHeatmap(viral.object, features = subacute.markers, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("subacute.markers_orig.ident_heatmap_counts.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = subacute.markers, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("subacute.markers_orig.ident_heatmap_logfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(viral.object) <- "Type"
feature.plot <- DotPlot(viral.object, features = inf.detected)
png(file=paste0("inf.detected-Type-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="Niche")
dev.off()
Idents(viral.object) <- "orig.ident"
library.averages.heatmap <- DoHeatmap(viral.object, features = inf.detected, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("inf.detected_orig.ident_heatmap_counts.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = inf.detected, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("inf.detected_orig.ident_heatmap_logfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(viral.object) <- "Type"
feature.plot <- DotPlot(viral.object, features = inf.notdetected)
png(file=paste0("inf.notdetected-Type-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="SARS-CoV-2 Transcripts") + scale_y_discrete(name ="Niche")
dev.off()
Idents(viral.object) <- "orig.ident"
library.averages.heatmap <- DoHeatmap(viral.object, features = inf.notdetected, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("inf.notdetected_orig.ident_heatmap_counts.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = inf.notdetected, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("inf.notdetected_orig.ident_heatmap_logfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")



Idents(viral.object) <- "Subtype"
unfiltered.count <- table(Idents(viral.object))
write.table(unfiltered.count, "Subtype.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(viral.object) <- "Type"
unfiltered.count <- table(Idents(viral.object))
write.table(unfiltered.count, "Type.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(viral.object) <- "Class"
unfiltered.count <- table(Idents(viral.object))
write.table(unfiltered.count, "Class.cellcounts.txt", sep="\t")
####################################################################################################
######################### Perform Pearson's Correlation Analysis of _all_ viral.only with _SARS.genes_   ########################################
########################################################################################################################
dir.create("Pearson_SARS.genes_all")
setwd("Pearson_SARS.genes_all")

Idents(seurat.object) <- "Type"

####Make a Pearson's Correlation Matrix script
## for "all" and "SARS.genes" from "viral.object" by "Type"#############
		##Mike Mariani UVM 2019
		##Hopefully useful for making a heatmap of Pearson correlation values
		##from the seurat scatterplots (where the title of the scatterplot)
		##is the Pearson's correlation value

library(Seurat)
library(corrplot)
library(Hmisc)
library(ggcorrplot)
library(data.table)

## What identity
Idents(viral.object) <- "Type"

## What identity is being subset?
# seurat.object3 <- subset(viral.object, subset = "all")
seurat.object3 <- viral.object
ncol(viral.object)


## What genes are were using for correlations?
features.input <- SARS.genes
features <- features.input


# set default assay from integrated back to Spatial. 
DefaultAssay(seurat.object3) <- "Spatial"
rows<-length(features)
cols<-length(features)    
x <- rep(NA, rows*cols)
cor.mat <- matrix(x, nrow=rows, ncol=cols)
scatter.list <- list()
count <- 1
for(i in 1:length(features))
{
  for(j in 1:length(features))     
  {
    scatter.list[[count]] <- Seurat::FeatureScatter(seurat.object3,
                                            as.character(features[i]),
                                            as.character(features[j]))
    cor.mat[i,j] <- scatter.list[[count]]$labels$title
    count <- count +1 
    print(paste0(features[i],",",features[j]))
  }
}

rownames(cor.mat) <- features
colnames(cor.mat) <- features
png(file=paste0("all_pearsons_heatmaps-hclus_all-SARS.genes.png"),
                res=300, 
                width=4500, 
                height=4500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 45,
         tl.cex=1)
dev.off()

png(file=paste0("all_pearsons_heatmaps-original_all-SARS.genes.png"),
                res=300, 
                width=4500, 
                height=4500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 45,
         tl.cex=1)
dev.off()

hmisc.cor.out <- Hmisc::rcorr(cor.mat)
hmisc.cor.out$P

##Compute a correlation matrix
corr <- round(cor(cor.mat), 1)
##Compute a matrix of correlation p-values
p.mat <- cor_pmat(cor.mat)
##Correlation matrix visualization
##Visualize the correlation matrix
##--------------------------------
##method = "square" (default)
png(file=paste0("all_pearsons_heatmaps_sig_hclustorder_all-SARS.genes.png"),
                res=300, 
                width=4500, 
                height=4500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         tl.cex=1)
dev.off()
##Correlation matrix visualization
png(file=paste0("all_pearsons_heatmaps_sig_origorder_all-SARS.genes.png"),
                res=300, 
                width=4500, 
                height=4500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         tl.cex=1)

dev.off()

setwd("..")
####################################################################################################
######################### Perform Pearson's Correlation Analysis of _all_ viral.only with _acute.markers_   ########################################
########################################################################################################################
dir.create("Pearson_acute.markers_all")
setwd("Pearson_acute.markers_all")
################################################################################################
####Make a Pearson's Correlation Matrix script
## for "all" and "acute.markers" from "viral.object" #############
		##Mike Mariani UVM 2019
		##Hopefully useful for making a heatmap of Pearson correlation values
		##from the seurat scatterplots (where the title of the scatterplot)
		##is the Pearson's correlation value

library(Seurat)
library(corrplot)
library(Hmisc)
library(ggcorrplot)
library(data.table)

## What identity
Idents(viral.object) <- "Class"

## What identity is being subset?
# seurat.object3 <- subset(viral.object, subset = "all")
seurat.object3 <- viral.object 
ncol(viral.object)


## What genes are were using for correlations?
features.input <- acute.markers
features <- features.input


# set default assay from integrated back to Spatial. 
DefaultAssay(seurat.object3) <- "Spatial"
rows<-length(features)
cols<-length(features)    
x <- rep(NA, rows*cols)
cor.mat <- matrix(x, nrow=rows, ncol=cols)
scatter.list <- list()
count <- 1
for(i in 1:length(features))
{
  for(j in 1:length(features))     
  {
    scatter.list[[count]] <- Seurat::FeatureScatter(seurat.object3,
                                            as.character(features[i]),
                                            as.character(features[j]))
    cor.mat[i,j] <- scatter.list[[count]]$labels$title
    count <- count +1 
    print(paste0(features[i],",",features[j]))
  }
}

rownames(cor.mat) <- features
colnames(cor.mat) <- features
png(file=paste0("all_pearsons_heatmaps-hclus_all-acute.markers.png"),
                res=300, 
                width=4500, 
                height=4500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 45,
         tl.cex=1)
dev.off()

png(file=paste0("all_pearsons_heatmaps-original_all-acute.markers.png"),
                res=300, 
                width=4500, 
                height=4500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 45,
         tl.cex=1)
dev.off()

hmisc.cor.out <- Hmisc::rcorr(cor.mat)
hmisc.cor.out$P

##Compute a correlation matrix
corr <- round(cor(cor.mat), 1)
##Compute a matrix of correlation p-values
p.mat <- cor_pmat(cor.mat)
##Correlation matrix visualization
##Visualize the correlation matrix
##--------------------------------
##method = "square" (default)
png(file=paste0("all_pearsons_heatmaps_sig_hclustorder_all-acute.markers.png"),
                res=300, 
                width=4500, 
                height=4500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         tl.cex=1)
dev.off()
##Correlation matrix visualization
png(file=paste0("all_pearsons_heatmaps_sig_origorder_all-acute.markers.png"),
                res=300, 
                width=4500, 
                height=4500)
corrplot(cor.mat, 
         method=c("square"),
         type = "full", 
         order = "original", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat=p.mat,
         sig.level = 0.05,
         insig="blank",
         tl.cex=1)

dev.off()
setwd("..")
################################################################################################



dir.create("SCT_var.genes_UMAP")
setwd("SCT_var.genes_UMAP")
viral.object <- SCTransform(viral.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes = FALSE, verbose = T)
DefaultAssay(viral.object) <- "SCT"
# viral.object <- NormalizeData(viral.object, normalization.method = "LogNormalize", scale.factor = 10000)
# viral.object <- ScaleData(viral.object, features = all.genes, assay="SCT")
# viral.object <- SCTransform(viral.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes=FALSE, verbose = T)
var.genes <- viral.object@assays[["SCT"]]@var.features
# viral.object<- ScaleData(viral.object, features = var.genes, assay = "SCT", verbose = TRUE)

## Dimensionality reduction, clustering, and visualization
viral.object <- RunPCA(viral.object, assay = "SCT", verbose = T)
viral.object <- FindNeighbors(viral.object, reduction = "pca", dims = 1:30)
viral.object <- FindClusters(viral.object, resolution = 0.6, verbose = T)
viral.object <- RunUMAP(viral.object, reduction = "pca", dims = 1:30)

p4 <- DimPlot(viral.object, reduction = "umap", group.by = "seurat_clusters", cols = mypal3, label=T, repel=T, raster=T) + labs(title = NULL, color='Cluster')
p4 

save.image("viral.only_clustered_v1.RData")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

####################################################################################################
## Subset viral transcriptomes of each sample
dir.create("orig.ident_viral.only_stats")
setwd("orig.ident_viral.only_stats")
SARS.genes <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")

viral.object <- DietSeurat(viral.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
DefaultAssay(viral.object) <- "Spatial"
viral.object.viral.genes <- DietSeurat(viral.object, assays="Spatial", counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs=NULL)
viral.object.viral.genes <- subset(viral.object.viral.genes, features = SARS.genes)


Idents(viral.object.viral.genes) <- "orig.ident"
viral.object.S16 <- subset(viral.object.viral.genes, idents = c("S16"))
ncol(viral.object.S16)
viral.object.S22 <- subset(viral.object.viral.genes, idents = c("S22"))
ncol(viral.object.S22)
viral.object.S23a <- subset(viral.object.viral.genes, idents = c("S23a"))
ncol(viral.object.S23a)
viral.object.S23b <- subset(viral.object.viral.genes, idents = c("S23b"))
ncol(viral.object.S23b)
	# 752 spots
# write.table(viral.object@assays[["Spatial"]]@counts, "viral.counts.txt", sep="\t")
# Save a table with cell counts per cluster

	## S16 1, S22 1, S23a 15, S23b 735
## Subset an object with only virus genes
write.table(viral.object.S16@assays[["Spatial"]]@counts, "viral.object.S16_viral.counts.txt", sep="\t")
write.table(viral.object.S22@assays[["Spatial"]]@counts, "viral.object.S22_viral.counts.txt", sep="\t")
write.table(viral.object.S23a@assays[["Spatial"]]@counts, "viral.object.S23a_viral.counts.txt", sep="\t")
write.table(viral.object.S23b@assays[["Spatial"]]@counts, "viral.object.S23b_viral.counts.txt", sep="\t")

write.table(viral.object.S16@assays[["Spatial"]]@data, "viral.object.S16_viral.counts_data.txt", sep="\t")
write.table(viral.object.S22@assays[["Spatial"]]@data, "viral.object.S22_viral.counts_data.txt", sep="\t")
write.table(viral.object.S23a@assays[["Spatial"]]@data, "viral.object.S23a_viral.counts_data.txt", sep="\t")
write.table(viral.object.S23b@assays[["Spatial"]]@data, "viral.object.S23b_viral.counts_data.txt", sep="\t")


Idents(viral.object.S16) <- "Type"
unfiltered.count <- table(Idents(viral.object.S16))
write.table(unfiltered.count, "viral.object.S16_Type_cellcounts.txt", sep="\t")

Idents(viral.object.S22) <- "Type"
unfiltered.count <- table(Idents(viral.object.S22))
write.table(unfiltered.count, "viral.object.S22_Type_cellcounts.txt", sep="\t")

Idents(viral.object.S23a) <- "Type"
unfiltered.count <- table(Idents(viral.object.S23a))
write.table(unfiltered.count, "viral.object.S23a_Type_cellcounts.txt", sep="\t")

Idents(viral.object.S23b) <- "Type"
unfiltered.count <- table(Idents(viral.object.S23b))
write.table(unfiltered.count, "viral.object.S23b_Type_cellcounts.txt", sep="\t")

setwd("..")
####################################################################################################


library(ggpubr)
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

p2 <- DimPlot(viral.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(viral.object, reduction = "umap", group.by = "Type",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
p2
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(viral.object, reduction = "umap", group.by = "Type",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Niche')
p2
ggsave("UMAP-Type-Niche.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(viral.object, reduction = "umap", group.by = "Subtype",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='Niche')
p2
ggsave("UMAP-Subtype-Niche.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")

p3 <- DimPlot(viral.object, reduction = "umap", group.by = "Phase", cols = mypal3, raster=FALSE) + labs(title = NULL, color='Cell Cycle')
ggsave("UMAP-Phase.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p3 <- DimPlot(viral.object, reduction = "umap", group.by = "orig.ident", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='SampleID')
ggsave("UMAP-orig.ident.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p4 <- DimPlot(viral.object, reduction = "umap", group.by = "Condition", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='Condition')
ggsave("UMAP-Condition.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")


p4 <- FeaturePlot(viral.object, reduction = "umap", features = c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10"), slot="counts", cols = c("lightgrey", "darkgreen", "lightgreen"), pt.size = 0.5) 
ggsave("UMAP-SARS.genes_rawcounts.pdf", plot = p4, device = "pdf", width = 11, height = 8, units = "in")
p4 <- FeaturePlot(viral.object, reduction = "umap", features = c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10"), cols = c("lightgrey", "darkblue", "lightblue"), pt.size = 0.5, keep.scale='feature') 
ggsave("UMAP-SARS.genes_logcounts.pdf", plot = p4, device = "pdf", width = 11, height = 8, units = "in")
p4 <- VlnPlot(viral.object, features = c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10"), group.by = 'seurat_clusters', assay = "Spatial", slot="counts", log=F)
ggsave("VlnPlot-SARS.genes_raw.pdf", plot = p4, device = "pdf", width = 11, height = 8, units = "in")
p4 <- VlnPlot(viral.object, features = c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10"), group.by = 'seurat_clusters')
ggsave("VlnPlot-SARS.genes_log.pdf", plot = p4, device = "pdf", width = 11, height = 8, units = "in")


p4 <- FeatureScatter(viral.object, feature1 = "ORF10", feature2 = "CXCL10")
ggsave("FeatureScatter-ORF10_CXCL10.pdf", plot = p4, device = "pdf", width = 8, height = 11, units = "in")
viral.object.c1 <-  subset(viral.object, identity = "1")

Idents(viral.object) <- "seurat_clusters"
cluster.averages <- AverageExpression(viral.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=SARS.genes)
write.csv(cluster.averages@assays[["Spatial"]]@counts, file = "cluster_viralmeancounts.csv")
######### Use this to determine the sum of mean counts of viral transcripts per cluster
### High to low sum of viral mean counts : c1>c4>c2>c0>c3

library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("Viral.only-Class.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("Viral.only-Class.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")

SARS.receptors <- c("ACE2", "TMPRSS2", "CTSB", "CTSL", "HBB", "HBA", "SRY", "XIST","ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
SARS.genes.receptors <- union(SARS.genes, SARS.receptors)
Idents(viral.object) <- "orig.ident"
library.averages.heatmap <- DoHeatmap(viral.object, features = SARS.genes.receptors, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("Ace2.SARS.genes.receptors-HB-Sex_orig.ident_heatmap_counts_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = SARS.genes.receptors, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("Ace2.SARS.genes.receptors-HB-Sex_orig.ident_heatmap_logfc_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(viral.object) <- "seurat_clusters"
library.averages.heatmap <- DoHeatmap(viral.object, features = SARS.genes.receptors, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("Ace2.SARS.genes.receptors-HB-Sex_seurat_clusters_heatmap_counts_seurat_clusters.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = SARS.genes.receptors, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("Ace2.SARS.genes.receptors-HB-Sex_seurat_clusters_heatmap_logfc_seurat_clusters.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

SARS.receptors <- c("ORF1AB-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10", "ACE2", "TMPRSS2", "CTSB", "CTSL", "HBB", "HBA")
Idents(viral.object) <- "seurat_clusters"
library.averages.heatmap <- DoHeatmap(viral.object, features = SARS.receptors, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors_seurat_clusters_heatmap_counts.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = SARS.receptors, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("SARS.genes.receptors_seurat_clusters_heatmap_logfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")




image.list <- c("S16.4")
p5 <- SpatialFeaturePlot(viral.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S16.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S17.5")
p5 <- SpatialFeaturePlot(viral.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S17.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S18.6")
p5 <- SpatialFeaturePlot(viral.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S18.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")

image.list <- c("S23a.11")
p5 <- SpatialFeaturePlot(viral.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S23a.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S23b.12")
p5 <- SpatialFeaturePlot(viral.object, features = SARS.genes, images=image.list, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S23b.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S24.13")

image.list <- c("S16.4")
p5 <- SpatialFeaturePlot(viral.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S16scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S17.5")
p5 <- SpatialFeaturePlot(viral.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S17scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S18.6")
p5 <- SpatialFeaturePlot(viral.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S18scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S22.10")
p5 <- SpatialFeaturePlot(viral.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S22scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S23a.11")
p5 <- SpatialFeaturePlot(viral.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S23ascaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S23b.12")
p5 <- SpatialFeaturePlot(viral.object, features = SARS.genes, images=image.list, min.cutoff=0, max.cutoff=9, ncol=5)
ggsave("viral.markers_SpatialFeaturePlots-cropped-S23bscaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
image.list <- c("S24.13")

Idents(viral.object) <- "seurat_clusters"
image.list <- c("S16.4", "S22.10", "S23a.11", "S23b.12")
plot2<- SpatialDimPlot(viral.object, images=image.list) # + patchwork::plot_layout(ncol = 4)
ggsave("SpatialFeaturePlot_seurat_clusters_cropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialDimPlot(viral.object, images=image.list, crop=FALSE) # + patchwork::plot_layout(ncol = 4)
ggsave("SpatialFeaturePlot_seurat_clusters_uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

Idents(viral.object) <- "Type"
image.list <- c("S16.4", "S22.10", "S23a.11", "S23b.12")
plot2<- SpatialDimPlot(viral.object, images=image.list) # + patchwork::plot_layout(ncol = 4)
ggsave("SpatialFeaturePlot_Type_cropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialDimPlot(viral.object, images=image.list, crop=FALSE) # + patchwork::plot_layout(ncol = 4)
ggsave("SpatialFeaturePlot_Type_uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")

Idents(viral.object) <- "Subtype"
image.list <- c("S16.4", "S22.10", "S23a.11", "S23b.12")
plot2<- SpatialDimPlot(viral.object, images=image.list) # + patchwork::plot_layout(ncol = 4)
ggsave("SpatialFeaturePlot_Subtype_cropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialDimPlot(viral.object, images=image.list, crop=FALSE) # + patchwork::plot_layout(ncol = 4)
ggsave("SpatialFeaturePlot_Subtype_uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")


Idents(viral.object) <- "Class"
image.list <- c("S16.4", "S22.10", "S23a.11", "S23b.12")
plot2<- SpatialDimPlot(viral.object, images=image.list) # + patchwork::plot_layout(ncol = 4)
ggsave("SpatialFeaturePlot_Class_cropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")
plot2<- SpatialDimPlot(viral.object, images=image.list, crop=FALSE) # + patchwork::plot_layout(ncol = 4)
ggsave("SpatialFeaturePlot_Class_uncropped.pdf", plot = plot2, device = "pdf", width = 8, height = 11, units = "in")



#examine UMAPs with qc metrics 
p5 <- FeaturePlot(viral.object, features = 'nFeature_Spatial')
p6 <- FeaturePlot(viral.object, features = 'nCount_Spatial')
p7 <- FeaturePlot(viral.object, features = 'percent.mt')
p8 <- FeaturePlot(viral.object, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots-Spatial.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")

##### Perform differential expression between seurat_clusters
dir.create("DE_seurat_clusters")
setwd("DE_seurat_clusters")
	## viral.object <- SCTransform(viral.object, method = "glmGamPoi", assay = "SCT", conserve.memory=TRUE, verbose = T)
# viral.object <- SCTransform(viral.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes = FALSE, verbose = T)
DefaultAssay(viral.object) <- "Spatial"
viral.object <- NormalizeData(viral.object, normalization.method = "LogNormalize", scale.factor = 10000)
viral.object <- ScaleData(viral.object, features = all.genes, assay="Spatial")

Idents(viral.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(viral.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_seurat_clusters_pos-0.693lnFC.txt", sep="\t")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("seurat_clusters_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("seurat_clusters_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("seurat_clusters_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

Idents(viral.object) <- "seurat_clusters"
library.averages.heatmap <- DoHeatmap(viral.object, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("clusters_heatmap_counts.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = unique.top2, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("clusters_heatmap_logfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

de_markers <- FindAllMarkers(viral.object, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_seurat_clusters_0.2lnFC.txt", sep="\t")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("seurat_clusters_top2.marker-DotPlot_0.2lnfc.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("seurat_clusters_top3.marker-DotPlot_0.2lnfc.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("seurat_clusters_top5.marker-DotPlot_0.2lnfc.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

Idents(viral.object) <- "seurat_clusters"
library.averages.heatmap <- DoHeatmap(viral.object, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("clusters_heatmap_counts_0.2lnfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = unique.top2, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("clusters_heatmap_logfc_0.2lnfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
## Slim down the viral.objects to save space. Don't need to keep the scale.data
# viral.object <- DietSeurat(viral.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

#Build a heirarchical clustering tree to see relationships between UMAP clusters
viral.object <- BuildClusterTree(viral.object, assay = "Spatial", features = SARS.genes) 
png(file=paste0("SARS.genes_clustertree_plot.png"),
                res=300, 
                width=1500, 
                height=2500)
PlotClusterTree(viral.object)
dev.off()


de_markers <- FindAllMarkers(viral.object, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0, logfc.threshold =  0)
write.table(de_markers, "DEGs_by_seurat_clusters_nofilter.txt", sep="\t")



setwd("..")


## Still not very good selective markers for each cluster.
	## probably bc these spatial transcriptomes are 1-5 cells per spot. 
	## look at immunology markers. 

##### Perform differential expression between Type
dir.create("DE_Type")
setwd("DE_Type")
	## viral.object <- SCTransform(viral.object, method = "glmGamPoi", assay = "SCT", conserve.memory=TRUE, verbose = T)
# viral.object <- SCTransform(viral.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes = FALSE, verbose = T)
DefaultAssay(viral.object) <- "Spatial"
viral.object <- NormalizeData(viral.object, normalization.method = "LogNormalize", scale.factor = 10000)
viral.object <- ScaleData(viral.object, features = all.genes, assay="Spatial")

Idents(viral.object) <- "Type"
de_markers <- FindAllMarkers(viral.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_Type_pos-0.693lnFC.txt", sep="\t")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("Type_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Niche")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("Type_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Niche")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("Type_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Niche")
dev.off()

Idents(viral.object) <- "Type"
library.averages.heatmap <- DoHeatmap(viral.object, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("clusters_heatmap_counts.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = unique.top2, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("clusters_heatmap_logfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

de_markers <- FindAllMarkers(viral.object, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_Type_0.2lnFC.txt", sep="\t")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("Type_top2.marker-DotPlot_0.2lnfc.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Niche")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("Type_top3.marker-DotPlot_0.2lnfc.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Niche")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("Type_top5.marker-DotPlot_0.2lnfc.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Niche")
dev.off()

Idents(viral.object) <- "Type"
library.averages.heatmap <- DoHeatmap(viral.object, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("clusters_heatmap_counts_0.2lnfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = unique.top2, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("clusters_heatmap_logfc_0.2lnfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
## Slim down the viral.objects to save space. Don't need to keep the scale.data
# viral.object <- DietSeurat(viral.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")


de_markers <- FindAllMarkers(viral.object, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0, logfc.threshold =  0)
write.table(de_markers, "DEGs_by_Type_nofilter.txt", sep="\t")

setwd("..")

##### Perform differential expression between Subtype
dir.create("DE_Subtype")
setwd("DE_Subtype")
	## viral.object <- SCTransform(viral.object, method = "glmGamPoi", assay = "SCT", conserve.memory=TRUE, verbose = T)
# viral.object <- SCTransform(viral.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes = FALSE, verbose = T)
DefaultAssay(viral.object) <- "Spatial"
viral.object <- NormalizeData(viral.object, normalization.method = "LogNormalize", scale.factor = 10000)
viral.object <- ScaleData(viral.object, features = all.genes, assay="Spatial")

Idents(viral.object) <- "Subtype"
de_markers <- FindAllMarkers(viral.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_Subtype_pos-0.693lnFC.txt", sep="\t")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("Subtype_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Niche")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("Subtype_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Niche")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("Subtype_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Niche")
dev.off()

Idents(viral.object) <- "Subtype"
library.averages.heatmap <- DoHeatmap(viral.object, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("clusters_heatmap_counts.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = unique.top2, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("clusters_heatmap_logfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

de_markers <- FindAllMarkers(viral.object, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_Subtype_0.2lnFC.txt", sep="\t")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("Subtype_top2.marker-DotPlot_0.2lnfc.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Niche")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("Subtype_top3.marker-DotPlot_0.2lnfc.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Niche")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(viral.object, features = unique.top2)
png(file=paste0("Subtype_top5.marker-DotPlot_0.2lnfc.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Niche")
dev.off()

Idents(viral.object) <- "Subtype"
library.averages.heatmap <- DoHeatmap(viral.object, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("clusters_heatmap_counts_0.2lnfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(viral.object, features = unique.top2, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("clusters_heatmap_logfc_0.2lnfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
## Slim down the viral.objects to save space. Don't need to keep the scale.data
# viral.object <- DietSeurat(viral.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")


de_markers <- FindAllMarkers(viral.object, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0, logfc.threshold =  0)
write.table(de_markers, "DEGs_by_Subtype_nofilter.txt", sep="\t")

setwd("..")


dir.create("immunology.clusters")
setwd("immunology.clusters")
Idents(viral.object) <- "seurat_clusters"
## use RNA for all DE analysis and plots
# viral.object <- ScaleData(viral.object, features = all.genes, assay = "SCT")
viral.object <- ScaleData(viral.object, features = all.genes, assay = "Spatial")
DefaultAssay(viral.object) <- "Spatial"


Idents(viral.object) <- "seurat_clusters"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("ITGAM", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_ITGAM_THBD.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(viral.object, "ITGAM")
p2 <- FeaturePlot(viral.object, features="ITGAM", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(viral.object, "THBD")
p2 <- FeaturePlot(viral.object, features="THBD", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_THBD.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(viral.object) <- "seurat_clusters"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("CD4", "IFNG"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD4_IFNG.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(viral.object, "CD4")
p2 <- FeaturePlot(viral.object, features="CD4", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD4.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(viral.object, "IFNG")
p2 <- FeaturePlot(viral.object, features="IFNG", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_IFNG.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(viral.object) <- "seurat_clusters"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("CD38", "PTPRC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD38_PTPRC.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(viral.object, "CD38")
p2 <- FeaturePlot(viral.object, features="CD38", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD38.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(viral.object, "PTPRC")
p2 <- FeaturePlot(viral.object, features="PTPRC", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_PTPRC.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(viral.object) <- "seurat_clusters"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD38_FPR2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(viral.object, "CD38")
p2 <- FeaturePlot(viral.object, features="CD38", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD38.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(viral.object, "FPR2")
p2 <- FeaturePlot(viral.object, features="FPR2", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_FPR2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")



## kernel density estimation 
# https://www.nature.com/articles/s41591-021-01329-2#Sec8
major.t <- c("CD4", "CD8A", "CCR7", "PTPRC", "SELL", "CD27", "CD38", "CD44", "CXCR5", "CD40LG", "CCR7", "FOXP3", "IKZF2")
gd.t<- c("TRGV9", "TRDV2")
MAIT.t <- c("TRAV10", "TRAJ18", "TRAV1-2")
NK.t <- c("CD3G", "FCGR3A", "NCAM1", "NCR1")
TH1.t <- c("IFNG", "TBX21", "TNF")
TH2.t <- c("GATA3", "IL4")
TH17.t <- c("RORC", "IL17A", "IL17F", "IL21")
## https://www.bioconductor.org/packages/release/bioc/vignettes/Nebulosa/inst/doc/nebulosa_seurat.html#1_Overview
  ## https://academic.oup.com/bioinformatics/article-abstract/37/16/2485/6103785
# BiocManager::install("Nebulosa")
library("Nebulosa")
Idents(viral.object) <- "seurat_clusters"

# SARS.genes <- c("ORF1AB-1-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
p1 <- plot_density(viral.object, "CD4")
p2 <- FeaturePlot(viral.object, "CD8A")
p5 <- CombinePlots(plots = list(p1, p2))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("UMAP_KernelDensity_CD4-CD8A.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots
p3 <- plot_density(viral.object, c("CD4", "CD8A"))
p5 <- p3 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_CD4-CD8A_Joint.pdf", plot = p5, device = "pdf", width = 5, height = 7, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("CD4", "CD8A"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_major.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("TRGV9", "TRDV2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_gd.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("TRAV10", "TRAJ18"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MAIT.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("CD3G", "FCGR3A","NCAM1", "NCR1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_NK.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")

	#	Error in linbin2D.ks(x, gpoints[[1]], gpoints[[2]], w = w) : 
	#	  NA/NaN/Inf in foreign function call (arg 10)
	#	In addition: Warning message:
	#	In ks.defaults(x = x, w = w, binned = binned, bgridsize = bgridsize,  :
	#	  Weights don't sum to sample size - they have been scaled accordingly


## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("IFNG", "TBX21", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH1.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 16, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("GATA3", "IL4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH2.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("RORC", "IL17A", "IL17F", "IL21"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH17.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")




## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(viral.object, features = lymphoid.lineage)
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
feature.plot <- DotPlot(viral.object, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(viral.object, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

macrophage.lineage <- c("CD14", "CD68", "CSF1R", "ADGRE1", "SIGLEC1", "CD38", "GPR18", "FPR2", "MYC", "EGR2")
feature.plot <- DotPlot(viral.object, features = macrophage.lineage)
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

feature.plot <- DotPlot(viral.object, features = macrophage.lineage.markers.2)
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

cluster.averages <- AverageExpression(viral.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=t.lineages)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, slot="counts", assay="Spatial")
ggsave("t.lineages_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, assay="Spatial")
ggsave("t.lineages_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(viral.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=macrophage.lineage.markers.2)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("macrophage.lineage.markers.2_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, assay="Spatial")
ggsave("macrophage.lineage.markers.2_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(viral.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=lymphoid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("lymphoid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, assay="Spatial")
ggsave("lymphoid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(viral.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=myeloid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("myeloid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, assay="Spatial")
ggsave("myeloid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
# viral.object <- DietSeurat(viral.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

############## Activation gating strategy https://www.nature.com/articles/s41590-021-01049-2/figures/10
## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("CD3D", "NCAM1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TorNK_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("CD4", "KLRB1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRB1-CD161.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("CD4", "CD38"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD38.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(viral.object, c("CD4", "CD69"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD69.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("CD4", "KLRK1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRK1-NKG2D.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("CD4", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("FOXP3", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-FOXP3-CD4.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("KIT", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("TNF", "TNFRSF8"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.TNFRSF8.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")


p4 <- plot_density(viral.object, c("CD14", "ITGAM"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("SIGLEC1", "FOLR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("SIGLEC1", "DDX3Y"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC-fetalsex_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
	## Placenta-associated maternal macrophages PAMMS (FOLR2 neg, HLA positive); HBCs FOLR2+HLA- ; Thomas et al. Naomi McGovern, JEM 2020
p4 <- plot_density(viral.object, c("FOLR2", "HLA-DRA"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMMs_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM1_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM2_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("ITGAM", "CD14"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M0.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("CD80", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("CD11B", "MYC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("CD163", "MRC1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2a.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("HLA-DRA", "IL6"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("CD14", "CD163"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2c.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("CD14", "CD1C"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(viral.object, c("CD14", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p1 <- plot_density(viral.object, "THBD")
p2 <- FeaturePlot(viral.object, "THBD")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(viral.object, "ITGB3")
p2 <- FeaturePlot(viral.object, "ITGB3")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_Megakaryocyte.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(viral.object, "CD1C")
p2 <- FeaturePlot(viral.object, "CD1C")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(viral.object, "IL6")
p2 <- FeaturePlot(viral.object, "IL6")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(viral.object, "TNF")
p2 <- FeaturePlot(viral.object, "TNF")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(viral.object, "MYC")
p2 <- FeaturePlot(viral.object, "MYC")
p5 <- CombinePlots(plots = list(p1, p2))


ggsave("UMAP_M2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(viral.object, "FOLR2")
p2 <- FeaturePlot(viral.object, "FOLR2")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_HBCs-FOLR2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(viral.object, "HLA-DRA")
p2 <- FeaturePlot(viral.object, "HLA-DRA")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_PAMM-HLA-DRA.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(viral.object, "CD33")
p2 <- FeaturePlot(viral.object, "CD33")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeliodProgenitor-M0.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p4 <- plot_density(viral.object, c("KIT", "CD33"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_MyeliodProgenitor-M0-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")



## http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_PI3K_AKT_MTOR_SIGNALING.html
inflammation <- list(c("ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "ADM", "ADORA2B", "ADRM1", "AHR", "APLNR", "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1", "BEST1", "BST2", "BTG2", "C3AR1", "C5AR1", "CALCRL", "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7", "CCRL2", "CD14", "CD40", "CD48", "CD55", "CD69", "CD70", "CD82", "CDKN1A", "CHST2", "CLEC5A", "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1", "CXCL10", "CXCL11", "CXCL6", "CXCL8", "CXCL9", "CXCR6", "CYBB", "DCBLD2", "EBI3", "EDN1", "EIF2AK2", "EMP3", "EREG", "F3", "FFAR2", "FPR1", "FZD5", "GABBR1", "GCH1", "GNA15", "GNAI3", "GP1BA", "GPC3", "GPR132", "GPR183", "HAS2", "HBEGF", "HIF1A", "HPN", "HRH1", "ICAM1", "ICAM4", "ICOSLG", "IFITM1", "IFNAR1", "IFNGR2", "IL10", "IL10RA", "IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5", "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2", "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV", "MEP1A", "MET", "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA", "NLRP3", "NMI", "NMUR1", "NOD2", "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2", "PCDH7", "PDE4B", "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2", "PTGER4", "PTGIR", "PTPRE", "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16", "RHOG", "RIPK2", "RNF144B", "ROS1", "RTP4", "SCARF1", "SCN1B", "SELE", "SELENOS", "SELL", "SEMA4D", "SERPINE1", "SGMS2", "SLAMF1", "SLC11A2", "SLC1A2", "SLC28A2", "SLC31A1", "SLC31A2", "SLC4A4", "SLC7A1", "SLC7A2", "SPHK1", "SRI", "STAB1", "TACR1", "TACR3", "TAPBP", "TIMP1", "TLR1", "TLR2", "TLR3", "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10", "TNFSF15", "TNFSF9", "TPBG", "VIP"))
viral.object <- AddModuleScore(object = viral.object, features = inflammation, ctrl = 5, name = 'inflammation.score')
p1 <- FeaturePlot(viral.object, features = 'inflammation.score1')
ggsave("FeaturePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'inflammation.score1')
ggsave("RidgePlot_inflam.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "inflammation.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "inflammation.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_inflammation.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

ifn.alpha <- list(c("ADAR", "B2M", "BATF2", "BST2", "C1S", "CASP1", "CASP8", "CCRL2", "CD47", "CD74", "CMPK2", "CMTR1", "CNP", "CSF1", "CXCL10", "CXCL11", "DDX60", "DHX58", "EIF2AK2", "ELF1", "EPSTI1", "GBP2", "GBP4", "GMPR", "HELZ2", "HERC6", "HLA-C", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3", "IL15", "IL4R", "IL7", "IRF1", "IRF2", "IRF7", "IRF9", "ISG15", "ISG20", "LAMP3", "LAP3", "LGALS3BP", "LPAR6", "LY6E", "MOV10", "MVB12A", "MX1", "NCOA7", "NMI", "NUB1", "OAS1", "OASL", "OGFR", "PARP12", "PARP14", "PARP9", "PLSCR1", "PNPT1", "PROCR", "PSMA3", "PSMB8", "PSMB9", "PSME1", "PSME2", "RIPK2", "RNF31", "RSAD2", "RTP4", "SAMD9", "SAMD9L", "SELL", "SLC25A28", "SP110", "STAT2", "TAP1", "TDRD7", "TENT5A", "TMEM140", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TRIM5", "TXNIP", "UBA7", "UBE2L6", "USP18", "WARS1"))
viral.object <- AddModuleScore(object = viral.object, features = ifn.alpha, ctrl = 5, name = 'ifn.alpha.score')
p1 <- FeaturePlot(viral.object, features = 'ifn.alpha.score1')
ggsave("FeaturePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'ifn.alpha.score1')
ggsave("RidgePlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "ifn.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "ifn.alpha.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_ifn.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

ifn.beta <- list(c("ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM", "BST2", "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86", "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1", "CMPK2", "CMTR1", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX58", "DDX60", "DHX58", "EIF2AK2", "EIF4E3", "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "GBP4", "GBP6", "GCH1", "GPR18", "GZMA", "HELZ2", "HERC6", "HIF1A", "HLA-A", "HLA-B", "HLA-DMA", "HLA-DQA1", "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IFNAR2", "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1", "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7", "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1", "METTL7B", "MT2A", "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12", "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP", "PNPT1", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1", "RIPK2", "RNF213", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP", "SERPING1", "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2", "SP110", "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4", "STAT1", "STAT2", "STAT3", "STAT4", "TAP1", "TAPBP", "TDRD7", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8", "VCAM1", "WARS1", "XAF1", "XCL1", "ZBP1", "ZNFX1"))
viral.object <- AddModuleScore(object = viral.object, features = ifn.beta, ctrl = 5, name = 'ifn.beta.score')
p1 <- FeaturePlot(viral.object, features = 'ifn.beta.score1')
ggsave("FeaturePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'ifn.beta.score1')
ggsave("RidgePlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "ifn.beta.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "ifn.beta.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_ifn.beta.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

apoptosis <- list(c("ADD1", "AIFM3", "ANKH", "ANXA1", "APP", "ATF3", "AVPR1A", "BAX", "BCAP31", "BCL10", "BCL2L1", "BCL2L10", "BCL2L11", "BCL2L2", "BGN", "BID", "BIK", "BIRC3", "BMF", "BMP2", "BNIP3L", "BRCA1", "BTG2", "BTG3", "CASP1", "CASP2", "CASP3", "CASP4", "CASP6", "CASP7", "CASP8", "CASP9", "CAV1", "CCNA1", "CCND1", "CCND2", "CD14", "CD2", "CD38", "CD44", "CD69", "CDC25B", "CDK2", "CDKN1A", "CDKN1B", "CFLAR", "CLU", "CREBBP", "CTH", "CTNNB1", "CYLD", "DAP", "DAP3", "DCN", "DDIT3", "DFFA", "DIABLO", "DNAJA1", "DNAJC3", "DNM1L", "DPYD", "EBP", "EGR3", "EMP1", "ENO2", "ERBB2", "ERBB3", "EREG", "ETF1", "F2", "F2R", "FAS", "FASLG", "FDXR", "FEZ1", "GADD45A", "GADD45B", "GCH1", "GNA15", "GPX1", "GPX3", "GPX4", "GSN", "GSR", "GSTM1", "GUCY2D", "H1-0", "HGF", "HMGB2", "HMOX1", "HSPB1", "IER3", "IFITM3", "IFNB1", "IFNGR1", "IGF2R", "IGFBP6", "IL18", "IL1A", "IL1B", "IL6", "IRF1", "ISG20", "JUN", "KRT18", "LEF1", "LGALS3", "LMNA", "LUM", "MADD", "MCL1", "MGMT", "MMP2", "NEDD9", "NEFH", "PAK1", "PDCD4", "PDGFRB", "PEA15", "PLAT", "PLCB2", "PLPPR4", "PMAIP1", "PPP2R5B", "PPP3R1", "PPT1", "PRF1", "PSEN1", "PSEN2", "PTK2", "RARA", "RELA", "RETSAT", "RHOB", "RHOT2", "RNASEL", "ROCK1", "SAT1", "SATB1", "SC5D", "SLC20A1", "SMAD7", "SOD1", "SOD2", "SPTAN1", "SQSTM1", "TAP1", "TGFB2", "TGFBR3", "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF12A", "TNFSF10", "TOP2A", "TSPO", "TXNIP", "VDAC2", "WEE1", "XIAP"))
viral.object <- AddModuleScore(object = viral.object, features = apoptosis, ctrl = 5, name = 'apoptosis.score')
p1 <- FeaturePlot(viral.object, features = 'apoptosis.score1')
ggsave("FeaturePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'apoptosis.score1')
ggsave("RidgePlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "apoptosis.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "apoptosis.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_apoptosis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tnf.alpha <- list(c("ABCA1", "ACKR3", "AREG", "ATF3", "ATP2B1", "B4GALT1", "B4GALT5", "BCL2A1", "BCL3", "BCL6", "BHLHE40", "BIRC2", "BIRC3", "BMP2", "BTG1", "BTG2", "BTG3", "CCL2", "CCL20", "CCL4", "CCL5", "CCN1", "CCND1", "CCNL1", "CCRL2", "CD44", "CD69", "CD80", "CD83", "CDKN1A", "CEBPB", "CEBPD", "CFLAR", "CLCF1", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL11", "CXCL2", "CXCL3", "CXCL6", "DDX58", "DENND5A", "DNAJB4", "DRAM1", "DUSP1", "DUSP2", "DUSP4", "DUSP5", "EDN1", "EFNA1", "EGR1", "EGR2", "EGR3", "EHD1", "EIF1", "ETS2", "F2RL1", "F3", "FJX1", "FOS", "FOSB", "FOSL1", "FOSL2", "FUT4", "G0S2", "GADD45A", "GADD45B", "GCH1", "GEM", "GFPT2", "GPR183", "HBEGF", "HES1", "ICAM1", "ICOSLG", "ID2", "IER2", "IER3", "IER5", "IFIH1", "IFIT2", "IFNGR2", "IL12B", "IL15RA", "IL18", "IL1A", "IL1B", "IL23A", "IL6", "IL6ST", "IL7R", "INHBA", "IRF1", "IRS2", "JAG1", "JUN", "JUNB", "KDM6B", "KLF10", "KLF2", "KLF4", "KLF6", "KLF9", "KYNU", "LAMB3", "LDLR", "LIF", "LITAF", "MAFF", "MAP2K3", "MAP3K8", "MARCKS", "MCL1", "MSC", "MXD1", "MYC", "NAMPT", "NFAT5", "NFE2L2", "NFIL3", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NINJ1", "NR4A1", "NR4A2", "NR4A3", "OLR1", "PANX1", "PDE4B", "PDLIM5", "PER1", "PFKFB3", "PHLDA1", "PHLDA2", "PLAU", "PLAUR", "PLEK", "PLK2", "PLPP3", "PMEPA1", "PNRC1", "PPP1R15A", "PTGER4", "PTGS2", "PTPRE", "PTX3", "RCAN1", "REL", "RELA", "RELB", "RHOB", "RIPK2", "RNF19B", "SAT1", "SDC4", "SERPINB2", "SERPINB8", "SERPINE1", "SGK1", "SIK1", "SLC16A6", "SLC2A3", "SLC2A6", "SMAD3", "SNN", "SOCS3", "SOD2", "SPHK1", "SPSB1", "SQSTM1", "STAT5A", "TANK", "TAP1", "TGIF1", "TIPARP", "TLR2", "TNC", "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFRSF9", "TNFSF9", "TNIP1", "TNIP2", "TRAF1", "TRIB1", "TRIP10", "TSC22D1", "TUBB2A", "VEGFA", "YRDC", "ZBTB10", "ZC3H12A", "ZFP36"))
viral.object <- AddModuleScore(object = viral.object, features = tnf.alpha, ctrl = 5, name = 'tnf.alpha.score')
p1 <- FeaturePlot(viral.object, features = 'tnf.alpha.score1')
ggsave("FeaturePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'tnf.alpha.score1')
ggsave("RidgePlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "tnf.alpha.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "tnf.alpha.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_tnf.alpha.png", plot = p1, device = "png", width = 7, height = 5, units = "in")


hypoxia.list <- c("ACKR3", "ADM", "ADORA2B", "AK4", "AKAP12", "ALDOA", "ALDOB", "ALDOC", "AMPD3", "ANGPTL4", "ANKZF1", "ANXA2", "ATF3", "ATP7A", "B3GALT6", "B4GALNT2", "BCAN", "BCL2", "BGN", "BHLHE40", "BNIP3L", "BRS3", "BTG1", "CA12", "CASP6", "CAV1", "CAVIN1", "CAVIN3", "CCN1", "CCN2", "CCN5", "CCNG2", "CDKN1A", "CDKN1B", "CDKN1C", "CHST2", "CHST3", "CITED2", "COL5A1", "CP", "CSRP2", "CXCR4", "DCN", "DDIT3", "DDIT4", "DPYSL4", "DTNA", "DUSP1", "EDN2", "EFNA1", "EFNA3", "EGFR", "ENO1", "ENO2", "ENO3", "ERO1A", "ERRFI1", "ETS1", "EXT1", "F3", "FAM162A", "FBP1", "FOS", "FOSL2", "FOXO3", "GAA", "GALK1", "GAPDH", "GAPDHS", "GBE1", "GCK", "GCNT2", "GLRX", "GPC1", "GPC3", "GPC4", "GPI", "GRHPR", "GYS1", "HAS1", "HDLBP", "HEXA", "HK1", "HK2", "HMOX1", "HOXB9", "HS3ST1", "HSPA5", "IDS", "IER3", "IGFBP1", "IGFBP3", "IL6", "ILVBL", "INHA", "IRS2", "ISG20", "JMJD6", "JUN", "KDELR3", "KDM3A", "KIF5A", "KLF6", "KLF7", "KLHL24", "LALBA", "LARGE1", "LDHA", "LDHC", "LOX", "LXN", "MAFF", "MAP3K1", "MIF", "MT1E", "MT2A", "MXI1", "MYH9", "NAGK", "NCAN", "NDRG1", "NDST1", "NDST2", "NEDD4L", "NFIL3", "NOCT", "NR3C1", "P4HA1", "P4HA2", "PAM", "PCK1", "PDGFB", "PDK1", "PDK3", "PFKFB3", "PFKL", "PFKP", "PGAM2", "PGF", "PGK1", "PGM1", "PGM2", "PHKG1", "PIM1", "PKLR", "PKP1", "PLAC8", "PLAUR", "PLIN2", "PNRC1", "PPARGC1A", "PPFIA4", "PPP1R15A", "PPP1R3C", "PRDX5", "PRKCA", "PYGM", "RBPJ", "RORA", "RRAGD", "S100A4", "SAP30", "SCARB1", "SDC2", "SDC3", "SDC4", "SELENBP1", "SERPINE1", "SIAH2", "SLC25A1", "SLC2A1", "SLC2A3", "SLC2A5", "SLC37A4", "SLC6A6", "SRPX", "STBD1", "STC1", "STC2", "SULT2B1", "TES", "TGFB3", "TGFBI", "TGM2", "TIPARP", "TKTL1", "TMEM45A", "TNFAIP3", "TPBG", "TPD52", "TPI1", "TPST2", "UGP2", "VEGFA", "VHL", "VLDLR", "WSB1", "XPNPEP1", "ZFP36", "ZNF292")
viral.object <- AddModuleScore(object = viral.object, features = hypoxia.list, ctrl = 5, name = 'hypoxia.score')
p1 <- FeaturePlot(viral.object, features = 'hypoxia.score1')
ggsave("FeaturePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'hypoxia.score1')
ggsave("RidgePlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "hypoxia.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "hypoxia.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_hypoxia.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

complement.list <- c("ACTN2", "ADAM9", "ADRA2B", "AKAP10", "ANG", "ANXA5", "APOA4", "APOBEC3F", "APOBEC3G", "APOC1", "ATOX1", "BRPF3", "C1QA", "C1QC", "C1R", "C1S", "C2", "C3", "C4BPB", "C9", "CA2", "CALM1", "CALM3", "CASP1", "CASP10", "CASP3", "CASP4", "CASP5", "CASP7", "CASP9", "CBLB", "CCL5", "CD36", "CD40LG", "CD46", "CD55", "CD59", "CDA", "CDH13", "CDK5R1", "CEBPB", "CFB", "CFH", "CLU", "COL4A2", "CP", "CPM", "CPQ", "CR1", "CR2", "CSRP1", "CTSB", "CTSC", "CTSD", "CTSH", "CTSL", "CTSO", "CTSS", "CTSV", "CXCL1", "DGKG", "DGKH", "DOCK10", "DOCK4", "DOCK9", "DPP4", "DUSP5", "DUSP6", "DYRK2", "EHD1", "ERAP2", "F10", "F2", "F3", "F5", "F7", "F8", "FCER1G", "FCN1", "FDX1", "FN1", "FYN", "GATA3", "GCA", "GMFB", "GNAI2", "GNAI3", "GNB2", "GNB4", "GNG2", "GNGT2", "GP1BA", "GP9", "GPD2", "GRB2", "GZMA", "GZMB", "GZMK", "HNF4A", "HPCAL4", "HSPA1A", "HSPA5", "IL6", "IRF1", "IRF2", "IRF7", "ITGAM", "ITIH1", "JAK2", "KCNIP2", "KCNIP3", "KIF2A", "KLK1", "KLKB1", "KYNU", "L3MBTL4", "LAMP2", "LAP3", "LCK", "LCP2", "LGALS3", "LGMN", "LIPA", "LRP1", "LTA4H", "LTF", "LYN", "MAFF", "ME1", "MMP12", "MMP13", "MMP14", "MMP15", "MMP8", "MSRB1", "MT3", "NOTCH4", "OLR1", "PCLO", "PCSK9", "PDGFB", "PDP1", "PFN1", "PHEX", "PIK3CA", "PIK3CG", "PIK3R5", "PIM1", "PLA2G4A", "PLA2G7", "PLAT", "PLAUR", "PLEK", "PLG", "PLSCR1", "PPP2CB", "PPP4C", "PRCP", "PRDM4", "PREP", "PRKCD", "PRSS3", "PRSS36", "PSEN1", "PSMB9", "RABIF", "RAF1", "RASGRP1", "RBSN", "RCE1", "RHOG", "RNF4", "S100A12", "S100A13", "S100A9", "SCG3", "SERPINA1", "SERPINB2", "SERPINC1", "SERPINE1", "SERPING1", "SH2B3", "SIRT6", "SPOCK2", "SRC", "STX4", "TFPI2", "TIMP1", "TIMP2", "TMPRSS6", "TNFAIP3", "USP14", "USP15", "USP16", "USP8", "VCPIP1", "WAS", "XPNPEP1", "ZEB1", "ZFPM2")
viral.object <- AddModuleScore(object = viral.object, features = complement.list, ctrl = 5, name = 'complement.score')
p1 <- FeaturePlot(viral.object, features = 'complement.score1')
ggsave("FeaturePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'complement.score1')
ggsave("RidgePlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "complement.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "complement.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_complement.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

angiogenesis.list <- c("APOH", "APP", "CCND2", "COL3A1", "COL5A2", "CXCL6", "FGFR1", "FSTL1", "ITGAV", "JAG1", "JAG2", "KCNJ8", "LPL", "LRPAP1", "LUM", "MSX1", "NRP1", "OLR1", "PDGFA", "PF4", "PGLYRP1", "POSTN", "PRG2", "PTK2", "S100A4", "SERPINA5", "SLCO2A1", "SPP1", "STC1", "THBD", "TIMP1", "TNFRSF21", "VAV2", "VCAN", "VEGFA", "VTN")
viral.object <- AddModuleScore(object = viral.object, features = angiogenesis.list, ctrl = 5, name = 'angiogenesis.score')
p1 <- FeaturePlot(viral.object, features = 'angiogenesis.score1')
ggsave("FeaturePlot_angiogenesis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'angiogenesis.score1')
ggsave("RidgePlot_angiogenesis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "angiogenesis.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_angiogenesis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "angiogenesis.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_angiogenesis.png", plot = p1, device = "png", width = 7, height = 5, units = "in")


il2 <- list(c("ABCB1", "ADAM19", "AGER", "AHCY", "AHNAK", "AHR", "ALCAM", "AMACR", "ANXA4", "APLP1", "ARL4A", "BATF", "BATF3", "BCL2", "BCL2L1", "BHLHE40", "BMP2", "BMPR2", "CA2", "CAPG", "CAPN3", "CASP3", "CCND2", "CCND3", "CCNE1", "CCR4", "CD44", "CD48", "CD79B", "CD81", "CD83", "CD86", "CDC42SE2", "CDC6", "CDCP1", "CDKN1C", "CISH", "CKAP4", "COCH", "COL6A1", "CSF1", "CSF2", "CST7", "CTLA4", "CTSZ", "CXCL10", "CYFIP1", "DCPS", "DENND5A", "DHRS3", "DRC1", "ECM1", "EEF1AKMT1", "EMP1", "ENO3", "ENPP1", "EOMES", "ETFBKMT", "ETV4", "F2RL2", "FAH", "FAM126B", "FGL2", "FLT3LG", "FURIN", "GABARAPL1", "GADD45B", "GALM", "GATA1", "GBP4", "GLIPR2", "GPR65", "GPR83", "GPX4", "GSTO1", "GUCY1B1", "HIPK2", "HK2", "HOPX", "HUWE1", "ICOS", "IFITM3", "IFNGR1", "IGF1R", "IGF2R", "IKZF2", "IKZF4", "IL10", "IL10RA", "IL13", "IL18R1", "IL1R2", "IL1RL1", "IL2RA", "IL2RB", "IL3RA", "IL4R", "IRF4", "IRF6", "IRF8", "ITGA6", "ITGAE", "ITGAV", "ITIH5", "KLF6", "LCLAT1", "LIF", "LRIG1", "LRRC8C", "LTB", "MAFF", "MAP3K8", "MAP6", "MAPKAPK2", "MUC1", "MXD1", "MYC", "MYO1C", "MYO1E", "NCOA3", "NCS1", "NDRG1", "NFIL3", "NFKBIZ", "NOP2", "NRP1", "NT5E", "ODC1", "P2RX4", "P4HA1", "PDCD2L", "PENK", "PHLDA1", "PHTF2", "PIM1", "PLAGL1", "PLEC", "PLIN2", "PLPP1", "PLSCR1", "PNP", "POU2F1", "PRAF2", "PRKCH", "PRNP", "PTCH1", "PTGER2", "PTH1R", "PTRH2", "PUS1", "RABGAP1L", "RGS16", "RHOB", "RHOH", "RNH1", "RORA", "RRAGD", "S100A1", "SCN9A", "SELL", "SELP", "SERPINB6", "SERPINC1", "SH3BGRL2", "SHE", "SLC1A5", "SLC29A2", "SLC2A3", "SLC39A8", "SMPDL3A", "SNX14", "SNX9", "SOCS1", "SOCS2", "SPP1", "SPRED2", "SPRY4", "ST3GAL4", "SWAP70", "SYNGR2", "SYT11", "TGM2", "TIAM1", "TLR7", "TNFRSF18", "TNFRSF1B", "TNFRSF21", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFSF10", "TNFSF11", "TRAF1", "TTC39B", "TWSG1", "UCK2", "UMPS", "WLS", "XBP1"))
viral.object <- AddModuleScore(object = viral.object, features = il2, ctrl = 5, name = 'il2.score')
p1 <- FeaturePlot(viral.object, features = 'il2.score1')
ggsave("FeaturePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'il2.score1')
ggsave("RidgePlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "il2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "il2.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_il2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

tgfb <- list(c("ACVR1", "APC", "ARID4B", "BCAR3", "BMP2", "BMPR1A", "BMPR2", "CDH1", "CDK9", "CDKN1C", "CTNNB1", "ENG", "FKBP1A", "FNTA", "FURIN", "HDAC1", "HIPK2", "ID1", "ID2", "ID3", "IFNGR2", "JUNB", "KLF10", "LEFTY2", "LTBP2", "MAP3K7", "NCOR2", "NOG", "PMEPA1", "PPM1A", "PPP1CA", "PPP1R15A", "RAB31", "RHOA", "SERPINE1", "SKI", "SKIL", "SLC20A1", "SMAD1", "SMAD3", "SMAD6", "SMAD7", "SMURF1", "SMURF2", "SPTBN1", "TGFB1", "TGFBR1", "TGIF1", "THBS1", "TJP1", "TRIM33", "UBE2D3", "WWTR1", "XIAP"))
viral.object <- AddModuleScore(object = viral.object, features = tgfb, ctrl = 5, name = 'tgfb.score')
p1 <- FeaturePlot(viral.object, features = 'tgfb.score1')
ggsave("FeaturePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'tgfb.score1')
ggsave("RidgePlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "tgfb.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "tgfb.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_tgfb.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

il6 <- list(c("A2M", "ACVR1B", "ACVRL1", "BAK1", "CBL", "CCL7", "CCR1", "CD14", "CD36", "CD38", "CD44", "CD9", "CNTFR", "CRLF2", "CSF1", "CSF2", "CSF2RA", "CSF2RB", "CSF3R", "CXCL1", "CXCL10", "CXCL11", "CXCL13", "CXCL3", "CXCL9", "DNTT", "EBI3", "FAS", "GRB2", "HAX1", "HMOX1", "IFNAR1", "IFNGR1", "IFNGR2", "IL10RB", "IL12RB1", "IL13RA1", "IL15RA", "IL17RA", "IL17RB", "IL18R1", "IL1B", "IL1R1", "IL1R2", "IL2RA", "IL2RG", "IL3RA", "IL4R", "IL6", "IL6ST", "IL7", "IL9R", "INHBE", "IRF1", "IRF9", "ITGA4", "ITGB3", "JUN", "LEPR", "LTB", "LTBR", "MAP3K8", "MYD88", "OSMR", "PDGFC", "PF4", "PIK3R5", "PIM1", "PLA2G2A", "PTPN1", "PTPN11", "PTPN2", "REG1A", "SOCS1", "SOCS3", "STAM2", "STAT1", "STAT2", "STAT3", "TGFB1", "TLR2", "TNF", "TNFRSF12A", "TNFRSF1A", "TNFRSF1B", "TNFRSF21", "TYK2"))
viral.object <- AddModuleScore(object = viral.object, features = il6, ctrl = 5, name = 'il6.score')
p1 <- FeaturePlot(viral.object, features = 'il6.score1')
ggsave("FeaturePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'il6.score1')
ggsave("RidgePlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "il6.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "il6.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_il6.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

# http://www.gsea-msigdb.org/gsea/msigdb/cards/COATES_MACROPHAGE_M1_VS_M2_UP.html
M1 <- c("ABCA9", "ABHD1", "ACSS1", "ACYP2", "ADA", "AIF1", "ALAD", "ARFGEF3", "ARSB", "ATP6V0E2", "AXL", "BCKDHB", "BLOC1S6", "CADM1", "CAP1", "CAPN5", "CBX6", "CD59", "CFH", "CLBA1", "CNRIP1", "COLEC12", "COMT", "CRIM1", "CXCL14", "CXCR4", "DST", "DYNLT1", "EMC1", "ENO2", "FAM124A", "FAM135A", "FAM9B", "FGD2", "FILIP1L", "GALNT11", "GATM", "GDA", "GJA1", "GLO1", "GNB4", "HAUS2", "HDDC3", "HLA-DQA1", "HMGN3", "KCNJ10", "LAMA3", "LCORL", "LYPLAL1", "MAF", "MALAT1", "MARCKSL1", "MARCO", "MSR1", "NAT8L", "NRCAM", "OCEL1", "OGFRL1", "P2RY13", "PIANP", "PIK3AP1", "PLAAT3", "PLBD1", "PLXDC2", "PPP2R5C", "PTGER3", "RAB10", "RAPSN", "RASAL2", "RCBTB2", "RCN1", "RFX3", "RPL14", "SFI1", "SLC35A1", "SLC7A7", "SLCO2B1", "SRD5A3", "TGFBI", "TIFAB", "TM7SF3", "TOR3A", "TTC3", "TUBB2B", "TXNIP", "ZNF727")
viral.object <- AddModuleScore(object = viral.object, features = M1, ctrl = 5, name = 'M1.score')
p1 <- FeaturePlot(viral.object, features = 'M1.score1')
ggsave("FeaturePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'M1.score1')
ggsave("RidgePlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "M1.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "M1.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M1.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

M2<- c("ADAM8", "AK4", "AKR1E2", "ALDOC", "ANG", "ANKRD37", "ANOS1", "ARG1", "ARSK", "ATG4D", "ATP6V0D2", "BNIP3", "BST1", "C1QB", "C1QC", "C1QTNF12", "C5orf34", "CBLB", "CD24", "CD300C", "CD300LD", "CD5L", "CHIA", "COL20A1", "CRIPT", "CTSK", "CTSL", "CYTIP", "EFHD2", "EGLN3", "ERO1A", "F10", "F7", "FAM177A1", "FAM199X", "FAM241A", "FBLIM1", "FLRT2", "GASK1B", "GDF15", "GPRC5B", "HILPDA", "HLA-A", "HLA-B", "HS6ST1", "HSPA1A", "HSPA1B", "IGF2R", "IGHM", "IL18BP", "ITGB3", "LBP", "LIN7C", "LRRC27", "MCOLN3", "MFGE8", "MGST2", "MYO1F", "PADI4", "PDXDC1", "PINK1", "PKDCC", "PRDX2", "PROCR", "PTGER2", "QRICH1", "RGCC", "RIMBP2", "RPGRIP1", "RRAS2", "SCD", "SH2B2", "SLC2A1", "SNHG6", "SOAT1", "THBS1", "TJP2", "TLR1", "TMEM267", "TULP4", "UCHL1", "VEGFA", "XDH")
viral.object <- AddModuleScore(object = viral.object, features = M2, ctrl = 5, name = 'M2.score')
p1 <- FeaturePlot(viral.object, features = 'M2.score1')
ggsave("FeaturePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'M2.score1')
ggsave("RidgePlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "M2.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "M2.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M2.png", plot = p1, device = "png", width = 7, height = 5, units = "in")


M2.list<- c("MYC", "EGR2", "CD163", "MRC1", "TGFB1", "VEGFA", "HLA-DRA", "IL6")
viral.object <- AddModuleScore(object = viral.object, features = M2.list, ctrl = 5, name = 'M2.list.score')
p1 <- FeaturePlot(viral.object, features = 'M2.list.score1')
ggsave("FeaturePlot_M2.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'M2.list.score1')
ggsave("RidgePlot_M2.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "M2.list.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M2.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "M2.list.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M2.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

M1.list<- c("CD38", "FPR2")
viral.object <- AddModuleScore(object = viral.object, features = M1.list, ctrl = 5, name = 'M1.list.score')
p1 <- FeaturePlot(viral.object, features = 'M1.list.score1')
ggsave("FeaturePlot_M1.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(viral.object, features = 'M1.list.score1')
ggsave("RidgePlot_M1.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(viral.object, features = "M1.list.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M1.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(viral.object, feature1 = "M1.list.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M1.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")


setwd("..")



dir.create("cell.counts")
setwd("cell.counts")

viral.object <- DietSeurat(viral.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

Idents(viral.object) <- "orig.ident"
unfiltered.count <- table(Idents(viral.object))
write.table(unfiltered.count, "sample.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(viral.object) <- "seurat_clusters"
unfiltered.count <- table(Idents(viral.object))
write.table(unfiltered.count, "cluster.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(viral.object) <- "Subtype"
unfiltered.count <- table(Idents(viral.object))
write.table(unfiltered.count, "Subtype.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(viral.object) <- "Type"
unfiltered.count <- table(Idents(viral.object))
write.table(unfiltered.count, "Type.cellcounts.txt", sep="\t")

Idents(viral.object) <- "Class"
unfiltered.count <- table(Idents(viral.object))
write.table(unfiltered.count, "Class.cellcounts.txt", sep="\t")



dir.create("seurat_clusters.Type")
setwd("seurat_clusters.Type")
## seurat_clusters.cluster into seurat_clusters.Type
Idents(viral.object)<-"seurat_clusters"
clusters.0.0 <- subset(viral.object, idents = c("0"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "0.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(viral.object)<-"seurat_clusters"
clusters.0.0 <- subset(viral.object, idents = c("1"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "1.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(viral.object)<-"seurat_clusters"
clusters.0.0 <- subset(viral.object, idents = c("2"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "2.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(viral.object)<-"seurat_clusters"
clusters.0.0 <- subset(viral.object, idents = c("3"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "3.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(viral.object)<-"seurat_clusters"
clusters.0.0 <- subset(viral.object, idents = c("4"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "4.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")


setwd("..")
# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("seurat_clusters.Type/"))
manifest<-manifest%>%mutate(name=gsub(".Type.cellcounts.txt","",value))
manifest
  
      ## These commands below are required for at least 1 of the data frames. Make sure to change the file pre-fix if necessary [AF_Epithelial]
df<-read.table("seurat_clusters.Type/0.Type.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type","0")
 df2<-df
df<-df%>%pivot_longer(cols = `0`,names_to = "seurat_clusters",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("seurat_clusters.Type/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Type",X)
  df<-df%>%pivot_longer(cols = X,names_to = "seurat_clusters",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of Type ## there are 3 seurat_clusterss x 6 Type subseurat_clusterss
  ## # levels(x=viral.object2) <- c("M0", "M0.M1", "M1", "M0.M2", "M2", "DC")
# df2%>%mutate(perc=percs_by_group(count,group = seurat_clusters))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])
# df7<-my_fxn(X = manifest$value[6])


final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)
#final<-full_join(final,df7)



final<-final%>%
  mutate(seurat_clusters=gsub(".Type.cellcounts.txt","",seurat_clusters))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = seurat_clusters))
#  pivot_wider(names_from = seurat_clusters,values_from = count)

write.table(final,"seurat_clusters.Type_CellCounts_merged.tsv",sep = "\t",row.names = F)

final
​
​
# df3<-sapply(manifest$value, my_fxn)
# 
# manifest<-manifest%>%
#   mutate(results=lapply(FUN = my_fxn(manifest$value)))
# 
#                             
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18

final$Type<-as.factor(final$Type)

final$seurat_clusters<-as.factor(final$seurat_clusters)
final
library(ggpubr)
ggbarplot(data = final,
          x = "seurat_clusters",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)

png(file=paste0("seurat_clusters.Type_diversity.plot.png"),
                res=300, 
                width=2000, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "seurat_clusters",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()




## 	https://dcl-wrangle.stanford.edu/pivot-advanced.html
final.0 <- select(final, c("Type", "seurat_clusters", "Percent"))
final.2<-final.0%>%pivot_wider(names_from=seurat_clusters, values_from=Percent)

write.table(final.2,"seurat_clusters.Type_pivot_wider.tsv",sep = "\t",row.names = F)

#########################################################################################
#########################################################################################
#########################################################################################

dir.create("Class.Type")
setwd("Class.Type")
## Class.cluster into Class.Type
Idents(viral.object)<-"Class"
levels(viral.object)
		# [1] "Placenta.Subacute" "Placenta.Acute" 
		#  new.metadata <- c("Control", "Control", "Control", "Respiratory.Asymptomatic", "Placenta.Subacute", "Respiratory.Asymptomatic", "Respiratory.Asymptomatic", "Respiratory.Asymptomatic", 
			# "Respiratory.Symptomatic", "Respiratory.Symptomatic", "Placenta.Subacute", "Placenta.Acute", "Placenta.Acute", "Control","Control","Control")


clusters.0.0 <- subset(viral.object, idents = c("Placenta.Subacute"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "Placenta.Subacute.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(viral.object)<-"Class"
clusters.0.0 <- subset(viral.object, idents = c("Placenta.Acute"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "Placenta.Acute.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")




setwd("..")
# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("Class.Type/"))
manifest<-manifest%>%mutate(name=gsub(".Type.cellcounts.txt","",value))
manifest
  
      ## These commands below are required for at least 1 of the data frames. Make sure to change the file pre-fix if necessary [AF_Epithelial]
df<-read.table("Class.Type/Placenta.Subacute.Type.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type","Placenta.Subacute")
 df2<-df
df<-df%>%pivot_longer(cols = `Placenta.Subacute`,names_to = "Class",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("Class.Type/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Type",X)
  df<-df%>%pivot_longer(cols = X,names_to = "Class",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of Type ## there are 3 Classs x 6 Type subClasss
  ## # levels(x=viral.object2) <- c("M0", "M0.M1", "M1", "M0.M2", "M2", "DC")
# df2%>%mutate(perc=percs_by_group(count,group = Class))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])

# df7<-my_fxn(X = manifest$value[6])


final<-full_join(df2,df3)

#final<-full_join(final,df7)



final<-final%>%
  mutate(Class=gsub(".Type.cellcounts.txt","",Class))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = Class))
#  pivot_wider(names_from = Class,values_from = count)

write.table(final,"Class.Type_CellCounts_merged.tsv",sep = "\t",row.names = F)

final
​
​
# df3<-sapply(manifest$value, my_fxn)
# 
# manifest<-manifest%>%
#   mutate(results=lapply(FUN = my_fxn(manifest$value)))
# 
#                             
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18

final$Type<-as.factor(final$Type)

final$Class<-as.factor(final$Class)
final
library(ggpubr)
ggbarplot(data = final,
          x = "Class",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)

png(file=paste0("Class.Type_diversity.plot.png"),
                res=300, 
                width=2000, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "Class",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()

## 	https://dcl-wrangle.stanford.edu/pivot-advanced.html
final.0 <- select(final, c("Type", "Class", "Percent"))
final.2<-final.0%>%pivot_wider(names_from=Class, values_from=Percent)

write.table(final.2,"Class.Type_pivot_wider.tsv",sep = "\t",row.names = F)

setwd("..")


#########################################################################################
#########################################################################################
#########################################################################################
## Do this ananlysis with seurat.object
setwd("/home/ebarrozo/visium/results/seurat_human_v2")
load("/home/ebarrozo/visium/results/seurat_human_v2/annotated/human_spatial_data-integrated-annotated_v1.RData")


setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1")
dir.create("cell.counts")
setwd("cell.counts")
## Check cell annotations
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", cols = mypal3, label=T, repel=T, raster=T) + labs(title = NULL, color='Niche')
p4   ## clusters are annotated
## Check seurat_clusters is still intact
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", cols = mypal3, label=T, repel=T, raster=T) + labs(title = NULL, color='Cluster')
p4   ## They are



## Add annotion data
new.metadata <- c("PV", "BP", "CAM", "PVBP", 
	"PVBP", "PVBP", "PVBP", "PVBP", 
	"PVBP", "PVBP", "PVBP", "PVBP", 
	"PVBP", "PVBP", "PVBP", "PVBP")
Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
	## [1] "S01"  "S03"  "S04"  "S15"  "S16"  "S17"  "S18"  "S19"  "S20"  "S21"  "S22"  "S23a" "S23b" "S24" "S25"  "S26" 
	## Any new annotation has to be in this order
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
## add cluster Subtypes as a metadata column
seurat.object$Site <- Idents(seurat.object)

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("Control", "Control", "Control", "SARS-CoV-2.Asymptomatic", "SARS-CoV-2.Symptomatic", "SARS-CoV-2.Asymptomatic", "SARS-CoV-2.Asymptomatic", "SARS-CoV-2.Asymptomatic", "SARS-CoV-2.Symptomatic", "SARS-CoV-2.Symptomatic", "SARS-CoV-2.Symptomatic", "SARS-CoV-2.Symptomatic", "SARS-CoV-2.Symptomatic", "Control","Control","Control")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Condition <- Idents(seurat.object)
Idents(seurat.object) <- "Condition"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("F", "F", "F", "F", 
	"M", "F", "M", "M", 
	"F", "F", "M", "F", 
	"F", "M","M","M")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$FetalSex <- Idents(seurat.object)
Idents(seurat.object) <- "FetalSex"

## New annotations 8/3/2022
Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
 new.metadata <- c("Control", "Control", "Control", "Respiratory.Asymptomatic", "Placenta.Subacute", "Respiratory.Asymptomatic", "Respiratory.Asymptomatic", "Respiratory.Asymptomatic", "Respiratory.Symptomatic", "Respiratory.Symptomatic", "Placenta.Subacute", "Placenta.Acute", "Placenta.Acute", "Control","Control","Control")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Class <- Idents(seurat.object)
Idents(seurat.object) <- "Class"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
 new.metadata <- c("Control", "Control", "Control", "NotDetected", "Detected", "Detected", "Detected", "NotDetected", "NotDetected", "NotDetected", "Detected", "Detected", "Detected", "Control","Control","Control")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Orthogonal.Detection <- Idents(seurat.object)
Idents(seurat.object) <- "Orthogonal.Detection"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
 new.metadata <- c("Control", "Control", "Control", "NotDetected", "Detected", "NotDetected", "NotDetected", "NotDetected", "NotDetected", "NotDetected", "Detected", "Detected", "Detected", "Control","Control","Control")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Visium.Detection <- Idents(seurat.object)
Idents(seurat.object) <- "Visium.Detection"

dir.create("Class.Type")
setwd("Class.Type")
## Class.cluster into Class.Type
Idents(seurat.object)<-"Class"
levels(seurat.object)
		# [1] "Placenta.Subacute" "Placenta.Acute" 
		#  new.metadata <- c("Control", "Control", "Control", "Respiratory.Asymptomatic", "Placenta.Subacute", "Respiratory.Asymptomatic", "Respiratory.Asymptomatic", "Respiratory.Asymptomatic", 
			# "Respiratory.Symptomatic", "Respiratory.Symptomatic", "Placenta.Subacute", "Placenta.Acute", "Placenta.Acute", "Control","Control","Control")


clusters.0.0 <- subset(seurat.object, idents = c("Placenta.Subacute"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "Placenta.Subacute.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"Class"
clusters.0.0 <- subset(seurat.object, idents = c("Placenta.Acute"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "Placenta.Acute.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"Class"
clusters.0.0 <- subset(seurat.object, idents = c("Control"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "Control.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"Class"
clusters.0.0 <- subset(seurat.object, idents = c("Respiratory.Asymptomatic"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "Respiratory.Asymptomatic.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"Class"
clusters.0.0 <- subset(seurat.object, idents = c("Respiratory.Symptomatic"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "Respiratory.Symptomatic.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

## These tables are located /home/ebarrozo/visium/results/seurat_human_v2/revisions_v1/Class.Type

setwd("..")
# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("Class.Type/"))
manifest<-manifest%>%mutate(name=gsub(".Type.cellcounts.txt","",value))
manifest
  
      ## These commands below are required for at least 1 of the data frames. Make sure to change the file pre-fix if necessary [AF_Epithelial]
df<-read.table("Class.Type/Placenta.Subacute.Type.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type","Placenta.Subacute")
 df2<-df
df<-df%>%pivot_longer(cols = `Placenta.Subacute`,names_to = "Class",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("Class.Type/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Type",X)
  df<-df%>%pivot_longer(cols = X,names_to = "Class",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of Type ## there are 3 Classs x 6 Type subClasss
  ## # levels(x=seurat.object2) <- c("M0", "M0.M1", "M1", "M0.M2", "M2", "DC")
# df2%>%mutate(perc=percs_by_group(count,group = Class))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])
# df7<-my_fxn(X = manifest$value[6])


final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)
#final<-full_join(final,df7)



final<-final%>%
  mutate(Class=gsub(".Type.cellcounts.txt","",Class))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = Class))
#  pivot_wider(names_from = Class,values_from = count)

## save data in the cell.counts folder
setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1/cell.counts")

write.table(final,"Class.Type_CellCounts_merged.tsv",sep = "\t",row.names = F)

final
​
​
# df3<-sapply(manifest$value, my_fxn)
# 
# manifest<-manifest%>%
#   mutate(results=lapply(FUN = my_fxn(manifest$value)))
# 
#                             
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18

final$Type<-as.factor(final$Type)

final$Class<-as.factor(final$Class)
final
library(ggpubr)
ggbarplot(data = final,
          x = "Class",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)

png(file=paste0("Class.Type_diversity.plot.png"),
                res=300, 
                width=2000, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "Class",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()

## 	https://dcl-wrangle.stanford.edu/pivot-advanced.html
final.0 <- select(final, c("Type", "Class", "Percent"))
final.2<-final.0%>%pivot_wider(names_from=Class, values_from=Percent)

write.table(final.2,"Class.Type_pivot_wider.tsv",sep = "\t",row.names = F)

setwd("..")

setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1/cell.counts")
dir.create("orig.ident.Type")
setwd("orig.ident.Type")
## orig.ident.cluster into orig.ident.Type
Idents(seurat.object)<-"orig.ident"
levels(seurat.object)
	# 	[1] "S01"  "S03"  "S04"  "S15"  "S16"  "S17"  "S18"  "S19"  "S20"  "S21"  "S22"  "S23a" "S23b" "S24" 
	#   [15] "S25"  "S26"

clusters.0.0 <- subset(seurat.object, idents = c("S01"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S01.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S03"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S03.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S04"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S04.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S15"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S15.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S16"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S16.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S17"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S17.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S18"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S18.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S19"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S19.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S20"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S20.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S21"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S21.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S22"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S22.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S23a"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S23a.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S23b"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S23b.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S24"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S24.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S25"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S25.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")

Idents(seurat.object)<-"orig.ident"
clusters.0.0 <- subset(seurat.object, idents = c("S26"))
Idents(clusters.0.0) <- "Type"
unfiltered.count <- table(Idents(clusters.0.0))
write.table(unfiltered.count, "S26.Type.cellcounts.txt", sep="\t")
rm("clusters.0.0")



setwd("..")
# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("orig.ident.Type/"))
manifest<-manifest%>%mutate(name=gsub(".Type.cellcounts.txt","",value))
manifest
  
      ## These commands below are required for at least 1 of the data frames. Make sure to change the file pre-fix if necessary [AF_Epithelial]
df<-read.table("orig.ident.Type/S01.Type.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type","S01")
 df2<-df
df<-df%>%pivot_longer(cols = `S01`,names_to = "orig.ident",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("orig.ident.Type/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Type",X)
  df<-df%>%pivot_longer(cols = X,names_to = "orig.ident",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of Type ## there are 3 orig.idents x 6 Type suborig.idents
  ## # levels(x=seurat.object2) <- c("M0", "M0.M1", "M1", "M0.M2", "M2", "DC")
# df2%>%mutate(perc=percs_by_group(count,group = orig.ident))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])
df7<-my_fxn(X = manifest$value[6])
df8<-my_fxn(X = manifest$value[7])
df9<-my_fxn(X = manifest$value[8])
df10<-my_fxn(X = manifest$value[9])
df11<-my_fxn(X = manifest$value[10])
df12<-my_fxn(X = manifest$value[11])
df13<-my_fxn(X = manifest$value[12])
df14<-my_fxn(X = manifest$value[13])
df15<-my_fxn(X = manifest$value[14])
df16<-my_fxn(X = manifest$value[15])
df17<-my_fxn(X = manifest$value[16])




final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)
final<-full_join(final,df7)
final<-full_join(final,df8)
final<-full_join(final,df9)
final<-full_join(final,df10)
final<-full_join(final,df11)
final<-full_join(final,df12)
final<-full_join(final,df13)
final<-full_join(final,df14)
final<-full_join(final,df15)
final<-full_join(final,df16)
final<-full_join(final,df17)




final<-final%>%
  mutate(orig.ident=gsub(".Type.cellcounts.txt","",orig.ident))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = orig.ident))
#  pivot_wider(names_from = orig.ident,values_from = count)

write.table(final,"orig.ident.Type_CellCounts_merged.tsv",sep = "\t",row.names = F)

final
​
​
# df3<-sapply(manifest$value, my_fxn)
# 
# manifest<-manifest%>%
#   mutate(results=lapply(FUN = my_fxn(manifest$value)))
# 
#                             
library(ggpubr)
# mypal<-get_palette("aaas",14)
    # Error: Insufficient values in manual scale. 18 needed but only 14 provided.
mypal<-get_palette("ucscgb",22)
    ## make sure you have enough colors, here 18

final$Type<-as.factor(final$Type)

final$orig.ident<-as.factor(final$orig.ident)
final
library(ggpubr)
ggbarplot(data = final,
          x = "orig.ident",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)

png(file=paste0("orig.ident.Type_diversity.plot.png"),
                res=300, 
                width=2000, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "orig.ident",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()

## 	https://dcl-wrangle.stanford.edu/pivot-advanced.html
final.0 <- select(final, c("Type", "orig.ident", "Percent"))
final.2<-final.0%>%pivot_wider(names_from=orig.ident, values_from=Percent)

write.table(final.2,"orig.ident.Type_pivot_wider.tsv",sep = "\t",row.names = F)

setwd("..")


setwd("..")


#########################################################################################
#########################################################################################
#########################################################################################




set.seed(seed=1)
setwd("/home/ebarrozo/visium/results")
library(dplyr)
library(Seurat)  # Seurat_4.0.3
library(RNAransform)
library(ggplot2)
library(cowplot)
library(SeuratDisk)
library(dplyr)
library(patchwork)
# BiocManager::install("glmGamPoi")
library(glmGamPoi)
library(ggpubr)


setwd("/home/ebarrozo/visium/results/seurat_human_v2/integration-v3/rcpa.integration")
load("human-rcpa-integrated-umap_v3.RData")
Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
 [1] "PE1"         "PE2"         "PE3"         "PE4"         "PN1"         "PN2"        
 [7] "PN3C"        "PN3P"        "PN4C"        "PN4P"        "SRR10166475" "SRR10166477"
[13] "SRR10166480" "SRR10166483" "SRR10166486" "SRR10166489" "SRR10166492" "SRR10166495"
[19] "SRR10166498" "SRR10166484" "SRR10166481" "SRR10166496" "SRR10166478" "SRR10166487"
[25] "SRR10166490" "SRR10166493" "SRR10166479" "SRR10166488" "SRR10166494" "SRR10166499"
[31] "SRR10166476" "SRR10166482" "SRR10166485" "SRR10166491" "SRR10166497" "SRR12192596"
[37] "SRR12192597" "SRR12192598" "SRR12192599" "C1"          "C2"          "C3"         
[43] "C4"          "P1"          "P2"          "G1"          "G2"          "SRR17065313"
[49] "SRR17065314" "SRR17065315" "SRR17065316" "SRR17065317" "SRR17065318" "SRR17065319"
[55] "SRR17065320" "SRR17065305" "SRR17065306" "SRR17065307" "SRR17065308" "SRR17065309"
[61] "SRR17065310" "SRR17065311" "SRR17065312" "S01"         "S03"         "S04"        
[67] "S15"         "S16"         "S17"         "S18"         "S19"         "S20"        
[73] "S21"         "S22"         "S23a"        "S23b"        "S24"         "S25"        
[79] "S26"    



list.1 <- levels(seurat.object)
setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1/atlas")
write.table(list.1, "atlas_orig.ident.txt", sep="\t")



new.metadata <- c("M", "M", "F", "M", "F", "F", "F", "F", "F", "M", "M", "M", "M", "M", "F", "F", "M", "F", "F", "F", "M", "M", "M", "M", "M", "M", "M", "M", "M", "F", "F", "F", "M", "F", "M", "F", "F", "M", "M", "F", "M", "M", "M", "M", "M", "F", "F", "F", "F", "M", "F", "M", "M", "F", "F", "M", "F", "M", "M", "M", "M", "M", "M", "F", "F", "M", "M", "F", "F", "M", "F", "F", "F", "F", "F", "F", "F", "F", "F")
write.table(new.metadata, "atlas_fetalsex.txt", sep="\t")

rm(seurat.object)



setwd("/home/ebarrozo/placenta-scrnaseq/phs001886.v3/results/seurat-v1")
load("wayne_data-integrated_umap-v1-annotated.RData")
    ## see phs001886.v3_seurat-analysis_v1.R for preprocessing
Idents(seurat.object) <- "SampleID"
levels(seurat.object)
 [1] "s1DB"              "s2DB"              "s3DB"              "s4DB"             
 [5] "s5DB"              "s6DB"              "s7DB"              "s8DB"             
 [9] "s9DB"              "CC7-1"             "CC7-2"             "CC8-1"            
[13] "CC8-2"             "s4P"               "s3P"               "s8P"              
[17] "s2P"               "s5P"               "s6P"               "s7P"              
[21] "C19C-PVBP-10"      "C19C-PVBP-12"      "C19C-PVBP-14"      "C19C-PVBP-16"     
[25] "C19C-PVBP-18"      "C19C-PVBP-20"      "C19C-PVBP-22"      "C19C-PVBP-24"     
[29] "HPL20289_R_PVBP_1" "HPL20289_R_PVBP_2" "s2W"               "s5W"              
[33] "s7W"               "s9W"               "s1W"               "s3W"              
[37] "s4W"               "s6W"               "s8W"               "C19C-CAM-09"      
[41] "C19C-CAM-11"       "C19C-CAM-13"       "C19C-CAM-15"       "C19C-CAM-17"      
[45] "C19C-CAM-19"       "C19C-CAM-21"       "C19C-CAM-23" 


c("s1DB"      ,        "s2DB"      ,        "s3DB"       ,       "s4DB" ,            
"s5DB"        ,      "s6DB"          ,    "s7DB"        ,      "s8DB"       ,      
"s9DB"        ,      "CC7-1"         ,    "CC7-2"        ,     "CC8-1"       ,     
"CC8-2"       ,      "s4P"           ,    "s3P"           ,    "s8P"          ,    
"s2P"          ,     "s5P"           ,    "s6P"            ,   "s7P"           ,   
"C19C-PVBP-10"  ,    "C19C-PVBP-12"  ,    "C19C-PVBP-14"    ,  "C19C-PVBP-16"   ,  
"C19C-PVBP-18"   ,   "C19C-PVBP-20"  ,    "C19C-PVBP-22"    ,  "C19C-PVBP-24"    , 
"HPL20289_R_PVBP_1", "HPL20289_R_PVBP_2", "s2W"            ,   "s5W"              ,
"s7W"           ,    "s9W"           ,    "s1W"            ,   "s3W"              ,
"s4W"           ,    "s6W"           ,    "s8W"            ,   "C19C-CAM-09"      ,
"C19C-CAM-11"   ,    "C19C-CAM-13"   ,    "C19C-CAM-15"    ,   "C19C-CAM-17"      ,
"C19C-CAM-19"   ,    "C19C-CAM-21"   ,    "C19C-CAM-23")

Idents(seurat.object) <- "Site"
levels(seurat.object)


Idents(seurat.object) <- "SampleID"
new.metadata <- c("DB"      ,        "DB"      ,        "DB"       ,       "DB" ,            
"DB"        ,      "DB"          ,    "DB"        ,      "DB"       ,      
"DB"        ,      "PVBP"         ,    "PVBP"        ,     "PVBP"       ,     
"PVBP"       ,      "PV"           ,    "PV"           ,    "PV"          ,    
"PV"          ,     "PV"           ,    "PV"            ,   "PV"           ,   
"PVBP"  ,    "PVBP"  ,    "PVBP"    ,  "PVBP"   ,  
"PVBP"   ,   "PVBP"  ,    "PVBP"    ,  "PVBP"    , 
"PVBP", "PVBP", "CAM"            ,   "CAM"              ,
"CAM"           ,    "CAM"           ,    "CAM"            ,   "CAM"              ,
"CAM"           ,    "CAM"           ,    "CAM"            ,   "CAM"      ,
"CAM"   ,    "CAM"   ,    "CAM"    ,   "CAM"      ,
"CAM"   ,    "CAM"   ,    "CAM")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Site <- Idents(seurat.object)

Idents(seurat.object) <- "SampleID"
new.metadata <- c("Control", "Control",  "Control", "Control",             
"Control", "Control",  "Control", "Control", 
"Control", "Control",  "Control", "Control", 
"Control", "Control",  "Control", "Control", 
"Control", "Control",  "Control", "Control",  
"SARS-CoV-2",    "SARS-CoV-2",    "SARS-CoV-2",  "SARS-CoV-2",  
"SARS-CoV-2",   "SARS-CoV-2",    "SARS-CoV-2",  "SARS-CoV-2", 
"Control", "Control",  "Control", "Control", 
"Control", "Control",  "Control", "Control", 
"Control", "Control",  "Control",   "SARS-CoV-2",
"SARS-CoV-2",    "SARS-CoV-2",    "SARS-CoV-2",   "SARS-CoV-2",
"SARS-CoV-2",    "SARS-CoV-2",   "SARS-CoV-2")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Condition <- Idents(seurat.object)

Idents(seurat.object) <- "SampleID"
new.metadata <- c("Pique-Regi-2019", "Pique-Regi-2019", "Pique-Regi-2019",  "Pique-Regi-2019",
"Pique-Regi-2019", "Pique-Regi-2019", "Pique-Regi-2019",  "Pique-Regi-2019",     
"Pique-Regi-2019",      "Pique-Regi-2020", "Pique-Regi-2020", "Pique-Regi-2020",      
"Pique-Regi-2020",      "Pique-Regi-2019", "Pique-Regi-2019",  "Pique-Regi-2019",    
"Pique-Regi-2019", "Pique-Regi-2019", "Pique-Regi-2019",  "Pique-Regi-2019",  
"Garcia-Flores-2022", "Garcia-Flores-2022", "Garcia-Flores-2022", "Garcia-Flores-2022",
"Garcia-Flores-2022", "Garcia-Flores-2022", "Garcia-Flores-2022", "Garcia-Flores-2022",
"Pique-Regi-2020", "Pique-Regi-2020", "Pique-Regi-2019",  "Pique-Regi-2019",
"Pique-Regi-2019", "Pique-Regi-2019", "Pique-Regi-2019",  "Pique-Regi-2019",
"Pique-Regi-2019", "Pique-Regi-2019",  "Pique-Regi-2019",   "Garcia-Flores-2022",
"Garcia-Flores-2022", "Garcia-Flores-2022", "Garcia-Flores-2022", "Garcia-Flores-2022",
"Garcia-Flores-2022", "Garcia-Flores-2022", "Garcia-Flores-2022")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$etal <- Idents(seurat.object)

Idents(seurat.object) <- "SampleID"
new.metadata <- c("scRNA-seq", "scRNA-seq", "scRNA-seq",  "scRNA-seq",
"scRNA-seq", "scRNA-seq", "scRNA-seq",  "scRNA-seq",     
"scRNA-seq",      "snRNA-seq", "snRNA-seq", "snRNA-seq",      
"snRNA-seq",      "scRNA-seq", "scRNA-seq",  "scRNA-seq",    
"scRNA-seq", "scRNA-seq", "scRNA-seq",  "scRNA-seq",  
"scRNA-seq", "scRNA-seq", "scRNA-seq", "scRNA-seq",
"scRNA-seq", "scRNA-seq", "scRNA-seq", "scRNA-seq",
"snRNA-seq", "snRNA-seq", "scRNA-seq",  "scRNA-seq",
"scRNA-seq", "scRNA-seq", "scRNA-seq",  "scRNA-seq",
"scRNA-seq", "scRNA-seq",  "scRNA-seq",   "scRNA-seq",
"scRNA-seq", "scRNA-seq", "scRNA-seq", "scRNA-seq",
"scRNA-seq", "scRNA-seq", "scRNA-seq")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Platform <- Idents(seurat.object)

## Remove the second trimester snRNA-seq data , "HPL20289_R_PVBP_1", "HPL20289_R_PVBP_2"
Idents(seurat.object) <- "SampleID"
ncol(seurat.object)
    # 209533

seurat.object <- subset(x = seurat.object, idents = c("s1DB", "s2DB", "s2W",  "s3DB",
"s4DB", "s4P",  "s5DB", "s3P", 
"s5W",  "s6DB", "s7DB", "s7W", 
"s8P",  "s9W",  "s1W",  "s2P", 
"s3W",  "s4W",  "s5P",  "s6P", 
"s6W",  "s7P",  "s8DB", "s8W", 
"s9DB", "C19C-CAM-09", "C19C-CAM-11", "C19C-CAM-13", 
"C19C-CAM-15", "C19C-CAM-17", "C19C-CAM-19", "C19C-CAM-21", "C19C-CAM-23", 
"C19C-PVBP-10", "C19C-PVBP-12", "C19C-PVBP-14", "C19C-PVBP-16", "C19C-PVBP-18", 
"C19C-PVBP-20", "C19C-PVBP-22", "C19C-PVBP-24", "CC7-1",     
"CC7-2", "CC8-1", "CC8-2"))
ncol(seurat.object)
    # 203848
levels(seurat.object)
 [1] "s1DB"         "s2DB"         "s3DB"         "s4DB"         "s5DB"         "s6DB"        
 [7] "s7DB"         "s8DB"         "s9DB"         "CC7-1"        "CC7-2"        "CC8-1"       
[13] "CC8-2"        "s4P"          "s3P"          "s8P"          "s2P"          "s5P"         
[19] "s6P"          "s7P"          "C19C-PVBP-10" "C19C-PVBP-12" "C19C-PVBP-14" "C19C-PVBP-16"
[25] "C19C-PVBP-18" "C19C-PVBP-20" "C19C-PVBP-22" "C19C-PVBP-24" "s2W"          "s5W"         
[31] "s7W"          "s9W"          "s1W"          "s3W"          "s4W"          "s6W"         
[37] "s8W"          "C19C-CAM-09"  "C19C-CAM-11"  "C19C-CAM-13"  "C19C-CAM-15"  "C19C-CAM-17" 
[43] "C19C-CAM-19"  "C19C-CAM-21"  "C19C-CAM-23" 

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
list.1
    ## Are they in the same order?
list.2 <- levels(seurat.object)
setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1/atlas")
write.table(list.2, "waynedata_orig.ident.txt", sep="\t")
    ## was not in order, so made atlas and fetal sex metadata into ascending order, then paired waynedata.orig.ident with orig.order, then made wayne data ascending, added ascending fetal sex data, then ordered wayne data back to orig.order for corrected fetal sex order
        ## The fix above was not correct so determined fetal sex manually below
Idents(seurat.object) <- "orig.ident"
library.averages.heatmap <- VlnPlot(seurat.object, features = "DDX3Y", assay="RNA", group.by="orig.ident")
library.averages.heatmap
ggsave("vlnplot_DDX3Y_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

Idents(seurat.object) <- "orig.ident"
new.metadata <- c("M", "M", "M", "M", "F", "F", "F", "F", "F", "M", "M", "M", "M", "M", "F", "F", "M", "F", "F", "F", "M", "M", "M", "M", "M", "M", "M", "M", "M", "F", "F", "M", "M", "M", "M", "F", "F", "M", "M", "F", "M", "M", "M", "M", "M")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$FetalSex <- Idents(seurat.object)


library.averages.heatmap <- VlnPlot(seurat.object, features = "DDX3Y", assay="RNA", group.by="FetalSex")
library.averages.heatmap
ggsave("vlnplot_DDX3Y_FetalSex.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
## confirmed correct

DefaultAssay(seurat.object) <- "RNA"
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs=NULL, assays=c("RNA"))

setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1/atlas")


Idents(seurat.object) <- "Condition"
levels(seurat.object)   # "Control", "Respiratory.Asymptomatic", "Placenta.Subacute", Placenta.Acute , Respiratory.Symptomatic"
seurat.object.x <- subset(x = seurat.object, idents = c("SARS-CoV-2"))
ncol(seurat.object.x)
    # 147906


Idents(seurat.object.x) <- "Platform"
levels(seurat.object.x)
		# "scRNA-seq" ## only scRNA-seq SARS-CoV-2 data, not snRNA-seq. 
seurat.object.y <- subset(x = seurat.object.x, idents = c("snRNA-seq"))
ncol(seurat.object.y)

Idents(seurat.object.x) <- "etal"
levels(seurat.object.x)
	# [1] "Garcia-Flores-2022"

setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1/atlas")

dir.create("DE_FetalSex_Waynedata-SARS")
setwd("DE_FetalSex_Waynedata-SARS")
Idents(seurat.object.x) <- "FetalSex"
levels(seurat.object.x)

Idents(seurat.object.x) <- "FetalSex"
levels(seurat.object.x)
unfiltered.count <- table(Idents(seurat.object.x))
unfiltered.count
write.table(unfiltered.count, "Wayne_fetalsex_cellcounts.txt", sep="\t")
#      M      F 
# 138436   9470  

seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes, assay="RNA")

dimorphism.list <- c("FCGRT", "FCGR1A", "FCGR3B","IFI6", "CXCL10", "OAS1", "CCL2", "MX1", "IL10", "IFNA1", "IFNA2" ,"IFNG", "CCL4", "CD163", "DDX3Y")
Idents(seurat.object) <- "FetalSex"
cluster.averages <- AverageExpression(seurat.object, assays="RNA", slot="counts", return.seurat = TRUE, features=dimorphism.list)
write.csv(cluster.averages@assays[["RNA"]]@counts, file = "dimorphism.list_meancounts.csv")
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = dimorphism.list, raster = FALSE, slot="counts", assay="RNA")
library.averages.heatmap
        ## DDX3Y upregulated in female column
        ## redo female/male assignments
    ## fetal sex fixed now
## looks good now
ggsave("dimorphism.list_mc_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = dimorphism.list, raster = FALSE, assay="RNA")
ggsave("dimorphism.list_FC_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = dimorphism.list, raster = FALSE, assay="RNA", slot="counts")
ggsave("dimorphism.list_counts_heatmap_fetalsex.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = dimorphism.list, raster = FALSE, assay="RNA")
ggsave("dimorphism.list_fold-change_heatmap_fetalsex.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")

de_markers.1 <- FindAllMarkers(seurat.object.x, features = dimorphism.list, assay = "RNA", only.pos = FALSE, min.pct = 0, logfc.threshold =  0)
write.table(de_markers.1, "dimorphism.list_DEGs_by_FetalSex_pos-nofilter.txt", sep="\t")

de_markers <- FindAllMarkers(seurat.object.x, features = all.genes, assay = "RNA", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_FetalSex_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object.x, assays="RNA", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
male.gene <- "DDX3Y"
top2 <- union(top2, male.gene)
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
png(file=paste0("FetalSex_top30.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
feature.plot

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
male.gene <- "DDX3Y"
top2 <- union(top2, male.gene)
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.x, features = unique.top2)
png(file=paste0("FetalSex_top50.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()

# seurat.object.x <- ScaleData(seurat.object.x, features = all.genes, assay="RNA")

top5.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = FALSE, angle=90, size=3)
top5.heatmap
ggsave("FetalSex_top50_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = FALSE, slot="counts", angle=90, size=3)
ggsave("FetalSex_top50_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = TRUE, slot="counts", angle=90, size=3)
ggsave("FetalSex_top50_markers_counts_heatmap_rasterTRUE.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.x, features = unique.top2, raster = TRUE, angle=90, size=3)
ggsave("FetalSex_top50_markers_logfc_heatmap_rasterTRUE.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()

## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
male.gene <- "DDX3Y"
top2 <- union(top2, male.gene)
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="RNA")
ggsave("FetalSex.Top50_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="RNA")
ggsave("FetalSex.Top50_FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")

#de_markers <- FindAllMarkers(seurat.object.x, features = all.genes, assay = "RNA", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
#write.table(de_markers, "DEGs_by_FetalSex_0.2lnFC.txt", sep="\t")
setwd("..")
#########################################################################################
#########################################################################################
## Determine the cell.type assocated with each marker, subset those cells, then see if M/F is DEG
dimorphism.list <- c("FCGRT", "FCGR1A", "FCGR3B","IFI6", "CXCL10", "OAS1", "CCL2", "MX1", "IL10", "IFNA1", "IFNA2" ,"IFNG", "CCL4", "CD163")

FCGRT - DC Macrophage
FCGR1A same
FCGR3B ND
IFI6 SYT VCT
CXCL10 DC 
OAS1 DC
CCL2 VEC PVC DC Macrophage DC2 Endothelial
IL10RB VCT  
IFNG NK CD8T
IFNA2 ND
CCL4 NK DC
CCL4 Macrophage 
CD163 Macrophage
#########################################################################################
seurat.object.x <- DietSeurat(seurat.object.x, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs=NULL)
Idents(seurat.object.x) <- "Type"
levels(seurat.object.x)
 [1] "CD8T"          "Fibroblast"    "NK"            "Macrophage"    "EVT"          
 [6] "VCT"           "DC2"           "Stromal"       "SmoothMuscle"  "SYT"          
[11] "B"             "Megakaryocyte"


seurat.object.y <- subset(x = seurat.object.x, idents = c( "DC2", "Macrophage", "NK", "CD8T", "VCT", "SYT"))

dir.create("DE_FetalSex_Waynedata-SARS_DCmacTypeSub")
setwd("DE_FetalSex_Waynedata-SARS_DCmacTypeSub")
levels(seurat.object.y)

seurat.object.y <- NormalizeData(seurat.object.y, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object.y <- ScaleData(seurat.object.y, features = all.genes, assay="RNA")

library.averages.heatmap <- DoHeatmap(seurat.object.y, features = dimorphism.list, raster = FALSE, assay="RNA", slot="counts")
ggsave("Type_dimorphism.list_counts_heatmap_fetalsex.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.y, features = dimorphism.list, raster = FALSE, assay="RNA")
ggsave("Type_dimorphism.list_fold-change_heatmap_fetalsex.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")

Idents(seurat.object.y) <- "Type"
cluster.averages <- AverageExpression(seurat.object.y, assays="RNA", slot="counts", return.seurat = TRUE, features=dimorphism.list)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = dimorphism.list, raster = FALSE, slot="counts", assay="RNA")
library.averages.heatmap
ggsave("Type_dimorphism.list_mc_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = dimorphism.list, raster = FALSE, assay="RNA")
ggsave("Type_dimorphism.list_FC_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
## Redo analysis above

Idents(seurat.object.y) <- "FetalSex"
levels(seurat.object.y)
unfiltered.count <- table(Idents(seurat.object.y))
unfiltered.count
write.table(unfiltered.count, "fetalsex_cellcounts.txt", sep="\t")
 #    M     F 
# 86176  2870 



de_markers.1 <- FindAllMarkers(seurat.object.y, features = dimorphism.list, assay = "RNA", only.pos = FALSE, min.pct = 0, logfc.threshold =  0)
write.table(de_markers.1, "dimorphism.list_DEGs_by_FetalSex_pos-nofilter.txt", sep="\t")

de_markers <- FindAllMarkers(seurat.object.y, features = all.genes, assay = "RNA", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_FetalSex_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object.y, assays="RNA", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
male.gene <- "DDX3Y"
top2 <- union(top2, male.gene)
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.y, features = unique.top2)
png(file=paste0("FetalSex_top30.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
feature.plot

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
male.gene <- "DDX3Y"
top2 <- union(top2, male.gene)
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.y, features = unique.top2)
png(file=paste0("FetalSex_top50.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()

# seurat.object.y <- ScaleData(seurat.object.y, features = all.genes, assay="RNA")

top5.heatmap <- DoHeatmap(seurat.object.y, features = unique.top2, raster = FALSE, angle=90, size=3)
top5.heatmap
ggsave("FetalSex_top50_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.y, features = unique.top2, raster = FALSE, slot="counts", angle=90, size=3)
ggsave("FetalSex_top50_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.y, features = unique.top2, raster = TRUE, slot="counts", angle=90, size=3)
ggsave("FetalSex_top50_markers_counts_heatmap_rasterTRUE.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.y, features = unique.top2, raster = TRUE, angle=90, size=3)
ggsave("FetalSex_top50_markers_logfc_heatmap_rasterTRUE.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()

## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
male.gene <- "DDX3Y"
top2 <- union(top2, male.gene)
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="RNA")
ggsave("FetalSex.Top50_Counts-MeanCounts_heatmap.pdf", plot = library.averages.heatmap, device = "pdf", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="RNA")
ggsave("FetalSex.Top50_FoldChange-MeanCounts_heatmap.pdf", plot = library.averages.heatmap, device = "pdf", width = 8, height = 10, units = "in")



dimorphism.list <- c("FCGRT", "FCGR1A", "FCGR3B","IFI6", "CXCL10", "OAS1", "CCL2", "MX1", "IL10", "IFNA1", "IFNA2" ,"IFNG", "CCL4", "CD163", "DDX3Y")
Idents(seurat.object) <- "FetalSex"
cluster.averages <- AverageExpression(seurat.object, assays="RNA", slot="counts", return.seurat = TRUE, features=dimorphism.list)
write.csv(cluster.averages@assays[["RNA"]]@counts, file = "dimorphism.list_meancounts.csv")
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = dimorphism.list, raster = FALSE, slot="counts", assay="RNA")
library.averages.heatmap
        ## DDX3Y upregulated in female column
        ## redo female/male assignments
    ## fetal sex fixed now
## looks good now
ggsave("dimorphism.list_mc_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = dimorphism.list, raster = FALSE, assay="RNA")
ggsave("dimorphism.list_FC_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.y, features = dimorphism.list, raster = FALSE, assay="RNA", slot="counts")
ggsave("dimorphism.list_counts_heatmap_fetalsex.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.y, features = dimorphism.list, raster = FALSE, assay="RNA")
ggsave("dimorphism.list_fold-change_heatmap_fetalsex.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")



#########################################################################################
#########################################################################################
#########################################################################################
## Slim down the seurat.object.xs to save space. Don't need to keep the scale.data
seurat.object.x <- DietSeurat(seurat.object.x, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

rm(seurat.object.x)
rm(seurat.object.detected)
#########################################################################################
#########################################################################################
#########################################################################################



setwd("/home/ebarrozo/visium/results/seurat_human_v2")
load("/home/ebarrozo/visium/results/seurat_human_v2/annotated/human_spatial_data-integrated-annotated_v1.RData")
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, assays=c("Spatial"))

## https://satijalab.org/seurat/articles/spatial_vignette.html
setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v1/atlas")
dir.create("crossval.SARS.sexual.dimorphism")
setwd("crossval.SARS.sexual.dimorphism")

## From venn diagram comparing Wayne SARS-CoV-2 data (not cell-type subset) vs Visium data (Detected samples subset)
    ## These were the sexually dimorphic DEGs cross-validated between single-cell and visium
crossval.SARS.sexual.dimorphism.genes <- c("SMIM3", "VKORC1", "PLP2", "TPI1", "ATP1B1", "RPS28", "RPL12", "CRYAB", "H3F3B", "RPS13", "MYL6", "NDUFA4L2", "BNIP3", "FAU", "ZFP36", "PGAM1", "LDHA", "RPS26", "CRLF1", "SELENOM", "PPDPF", "RPL39", "CD9", "TIMP1", "PDPN", "IFITM3", "SEC61G", "EMP3", "PID1", "EIF4EBP1", "S100A10", "ACTB", "FAM162A", "SH3BGRL3", "CCN1", "RPL41", "PDLIM3", "RPL11", "RPL8", "LMCD1", "MIF", "ACTA2", "F3", "GAPDH", "MT2A", "PGK1", "SNRPD2", "HSD17B1", "GDF15", "S100P", "C1QC", "CLDN4", "A2M", "FBN2", "CGA", "NOTUM", "MFAP5", "NEAT1", "SELENOP", "FBLN1", "CPVL", "AOC1", "PAPPA2", "VSIG4", "TINAGL1", "TFRC", "RPS4Y1", "FLT1", "DAB2", "EGFR", "STAB1", "PRG2", "ISM2", "SLC43A2", "GRN", "DDX17", "RNASE1", "GPX3", "KCTD12", "EPAS1", "SRSF5", "CSF1R", "SPINT2", "KRT23")

## These genes were the ones included in the significant KEGG pathways
KEGG.sexual.dimorphism.genes <- c("VSIG4", "A2M", "C1QC", "RPS26", "RPL41", "RPS28", "RPL12", "RPL11", "FAU", "RPL8", "RPL39", "RPS13", "LDHA", "TPI1", "PGAM1", "PGK1", "GAPDH","EIF4EBP1", "TIMP1")

image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
image.list <- c("S01")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S01.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S03.1")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S03.1.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S04.2")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S04.2.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S15.3")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S15.3.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")

image.list <- c("S16.4")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S16.4.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S17.5")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S17.5.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S18.6")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S18.6.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S19.7")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S19.7.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S20.8")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S20.8.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S21.9")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S21.9.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")

image.list <- c("S21.9")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S21.9.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S23a.11")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S23a.11.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S23b.12")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S23b.12.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S24.13")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S24.13.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S25.14")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S25.14.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S26.15")
p1<- SpatialPlot(seurat.object, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S26.15.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

seurat.object <- ScaleData(seurat.object, features = all.genes, assay="Spatial")
Idents(seurat.object) <- "FetalSex"
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=KEGG.sexual.dimorphism.genes)
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, slot="counts", assay="Spatial")
library.averages.heatmap
ggsave("KEGG.sexual.dimorphism.genes_mc_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_FC_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial", slot="counts")
ggsave("KEGG.sexual.dimorphism.genes_counts_heatmap_fetalsex.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_fold-change_heatmap_fetalsex.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")

Idents(seurat.object) <- "Type"
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=KEGG.sexual.dimorphism.genes)
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, slot="counts", assay="Spatial")
library.averages.heatmap
ggsave("KEGG.sexual.dimorphism.genes_mc_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_FC_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial", slot="counts")
ggsave("KEGG.sexual.dimorphism.genes_counts_heatmap_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_fold-change_heatmap_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")

Idents(seurat.object) <- "Condition"
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=KEGG.sexual.dimorphism.genes)
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, slot="counts", assay="Spatial")
library.averages.heatmap
ggsave("KEGG.sexual.dimorphism.genes_mc_heatmap_Condition.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_FC_heatmap_Condition.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial", slot="counts")
ggsave("KEGG.sexual.dimorphism.genes_counts_heatmap_Condition.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_fold-change_heatmap_Condition.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")


Idents(seurat.object) <- "Visium.Detection"
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=KEGG.sexual.dimorphism.genes)
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, slot="counts", assay="Spatial")
library.averages.heatmap
ggsave("KEGG.sexual.dimorphism.genes_mc_heatmap_Visium.Detection.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_FC_heatmap_Visium.Detection.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial", slot="counts")
ggsave("KEGG.sexual.dimorphism.genes_counts_heatmap_Visium.Detection.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_fold-change_heatmap_Visium.Detection.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")


Idents(seurat.object) <- "Orthogonal.Detection"
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=KEGG.sexual.dimorphism.genes)
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, slot="counts", assay="Spatial")
library.averages.heatmap
ggsave("KEGG.sexual.dimorphism.genes_mc_heatmap_Orthogonal.Detection.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_FC_heatmap_Orthogonal.Detection.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial", slot="counts")
ggsave("KEGG.sexual.dimorphism.genes_counts_heatmap_Orthogonal.Detection.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_fold-change_heatmap_Orthogonal.Detection.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")



Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("Control", "Control", "Control", "SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-2", "Control","Control","Control")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Status <- Idents(seurat.object)
Idents(seurat.object) <- "Status"

Idents(seurat.object) <- "Status"
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=KEGG.sexual.dimorphism.genes)
write.csv(cluster.averages@assays[["Spatial"]]@counts, file = "KEGG.sexual.dimorphism.genes_meancounts.csv")
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, slot="counts", assay="Spatial")
library.averages.heatmap
ggsave("KEGG.sexual.dimorphism.genes_mc_heatmap_Status.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_FC_heatmap_Status.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial", slot="counts")
ggsave("KEGG.sexual.dimorphism.genes_counts_heatmap_Status.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_fold-change_heatmap_Status.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")



dir.create("DE_Status")
setwd("DE_Status")

Idents(seurat.object) <- "Status"
levels(seurat.object)
unfiltered.count <- table(Idents(seurat.object))
unfiltered.count
write.table(unfiltered.count, "Status_cellcounts.txt", sep="\t")
 #    M     F 
# 86176  2870 


dimorphism.list <- c("FCGRT", "FCGR1A", "FCGR3B","IFI6", "CXCL10", "OAS1", "CCL2", "MX1", "IL10", "IFNA1", "IFNA2" ,"IFNG", "CCL4", "CD163")
Idents(seurat.object) <- "Status"

de_markers.1 <- FindAllMarkers(seurat.object, features = KEGG.sexual.dimorphism.genes, assay = "Spatial", only.pos = FALSE, min.pct = 0, logfc.threshold =  0)
write.table(de_markers.1, "KEGG.sexual.dimorphism.genes_DEGs_by_Status_pos-nofilter.txt", sep="\t")
de_markers.1 <- FindAllMarkers(seurat.object, features = dimorphism.list, assay = "Spatial", only.pos = FALSE, min.pct = 0, logfc.threshold =  0)
write.table(de_markers.1, "dimorphism.list_DEGs_by_Status_pos-nofilter.txt", sep="\t")

de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_Status_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
male.gene <- "DDX3Y"
top2 <- union(top2, male.gene)
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Status_top30.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Status")
dev.off()
feature.plot

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
male.gene <- "DDX3Y"
top2 <- union(top2, male.gene)
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Status_top50.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Status")
dev.off()

# seurat.object <- ScaleData(seurat.object, features = all.genes, assay="Spatial")

top5.heatmap <- DoHeatmap(seurat.object, features = unique.top2, raster = FALSE, angle=90, size=3)
top5.heatmap
ggsave("Status_top50_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object, features = unique.top2, raster = FALSE, slot="counts", angle=90, size=3)
ggsave("Status_top50_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object, features = unique.top2, raster = TRUE, slot="counts", angle=90, size=3)
ggsave("Status_top50_markers_counts_heatmap_rasterTRUE.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object, features = unique.top2, raster = TRUE, angle=90, size=3)
ggsave("Status_top50_markers_logfc_heatmap_rasterTRUE.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()

Idents(seurat.object) <- "Type"
cluster.averages <- AverageExpression(seurat.object, assays="Spatial", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
library.averages.heatmap
ggsave("Type_unique.top2_mc_heatmap_Status.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("Type_unique.top2_FC_heatmap_Status.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
male.gene <- "DDX3Y"
top2 <- union(top2, male.gene)
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("Status.Top50_Counts-MeanCounts_heatmap.pdf", plot = library.averages.heatmap, device = "pdf", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("Status.Top50_FoldChange-MeanCounts_heatmap.pdf", plot = library.averages.heatmap, device = "pdf", width = 8, height = 10, units = "in")


setwd("..")


dir.create("DE_Detected_KEGG.Crossval")
setwd("DE_Detected_KEGG.Crossval")


Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
 new.metadata <- c("Control", "Control", "Control", "NotDetected", "Detected", "NotDetected", "NotDetected", "NotDetected", "NotDetected", "NotDetected", "Detected", "Detected", "Detected", "Control","Control","Control")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Visium.Detection <- Idents(seurat.object)
Idents(seurat.object) <- "Visium.Detection"


Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("F", "F", "F", "F", 
	"M", "F", "M", "M", 
	"F", "F", "M", "F", 
	"F", "M","M","M")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$FetalSex <- Idents(seurat.object)
Idents(seurat.object) <- "FetalSex"


Idents(seurat.object) <- "Visium.Detection"
seurat.object.detected <- subset(x = seurat.object, idents = c("Detected"))



Idents(seurat.object.detected) <- "FetalSex"
levels(seurat.object.detected)
unfiltered.count <- table(Idents(seurat.object.detected))
unfiltered.count
write.table(unfiltered.count, "FetalSex_cellcounts.txt", sep="\t")
 #    F    M 
	# 1211 2173 

## v1 included log2 fold-change >0.693 not ln fold-change >0.693. 
KEGG.sexual.dimorphism.genes.v2 <- c("RPS4Y1", "RPS13", "EGFR", "TIMP1", "LDHA", "GAPDH", "TPI1", "RPL39", "PGK1", "RPS13", "RPL41")
KEGG.sexual.dimorphism.genes <- KEGG.sexual.dimorphism.genes.v2

dimorphism.list <- c("FCGRT", "FCGR1A", "FCGR3B","IFI6", "CXCL10", "OAS1", "CCL2", "MX1", "IL10", "IFNA1", "IFNA2" ,"IFNG", "CCL4", "CD163")
Idents(seurat.object.detected) <- "FetalSex"

de_markers.1 <- FindAllMarkers(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, assay = "Spatial", only.pos = FALSE, min.pct = 0, logfc.threshold =  0)
write.table(de_markers.1, "KEGG.sexual.dimorphism.genes_DEGs_by_FetalSex_pos-nofilter.txt", sep="\t")
de_markers.1 <- FindAllMarkers(seurat.object.detected, features = dimorphism.list, assay = "Spatial", only.pos = FALSE, min.pct = 0, logfc.threshold =  0)
write.table(de_markers.1, "dimorphism.list_DEGs_by_FetalSex_pos-nofilter.txt", sep="\t")


# seurat.object.detected <- ScaleData(seurat.object.detected, features = all.genes, assay="Spatial")
Idents(seurat.object.detected) <- "FetalSex"
cluster.averages <- AverageExpression(seurat.object.detected, assays="Spatial", slot="counts", return.seurat = TRUE, features=KEGG.sexual.dimorphism.genes)
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, slot="counts", assay="Spatial")
library.averages.heatmap
ggsave("KEGG.sexual.dimorphism.genes_mc_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_FC_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial", slot="counts")
ggsave("KEGG.sexual.dimorphism.genes_counts_heatmap_fetalsex.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_fold-change_heatmap_fetalsex.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")

Idents(seurat.object.detected) <- "Type"
cluster.averages <- AverageExpression(seurat.object.detected, assays="Spatial", slot="counts", return.seurat = TRUE, features=KEGG.sexual.dimorphism.genes)
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, slot="counts", assay="Spatial")
library.averages.heatmap
ggsave("KEGG.sexual.dimorphism.genes_mc_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_FC_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial", slot="counts")
ggsave("KEGG.sexual.dimorphism.genes_counts_heatmap_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, raster = FALSE, assay="Spatial")
ggsave("KEGG.sexual.dimorphism.genes_fold-change_heatmap_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")

setwd("..")

## 	REDO SPATIAL PLOTS ABOVE
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
image.list <- c("S01")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S01.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S03.1")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S03.1.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S04.2")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S04.2.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S15.3")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S15.3.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")

image.list <- c("S16.4")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S16.4.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S17.5")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S17.5.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S18.6")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S18.6.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S19.7")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S19.7.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S20.8")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S20.8.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S21.9")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S21.9.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")

image.list <- c("S21.9")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S21.9.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S23a.11")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S23a.11.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S23b.12")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S23b.12.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S24.13")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S24.13.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S25.14")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S25.14.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S26.15")
p1<- SpatialPlot(seurat.object.detected, features = KEGG.sexual.dimorphism.genes, images = image.list) + patchwork::plot_layout(ncol = 4)
p1
ggsave("KEGG.sexual.dimorphism.genes_S26.15.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")





de_markers <- FindAllMarkers(seurat.object.detected, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_FetalSex_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object.detected, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
male.gene <- "DDX3Y"
top2 <- union(top2, male.gene)
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.detected, features = unique.top2)
png(file=paste0("FetalSex_top30.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
feature.plot

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
male.gene <- "DDX3Y"
top2 <- union(top2, male.gene)
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.detected, features = unique.top2)
png(file=paste0("FetalSex_top50.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()

# seurat.object.detected <- ScaleData(seurat.object.detected, features = all.genes, assay="Spatial")

top5.heatmap <- DoHeatmap(seurat.object.detected, features = unique.top2, raster = FALSE, angle=90, size=3)
top5.heatmap
ggsave("FetalSex_top50_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.detected, features = unique.top2, raster = FALSE, slot="counts", angle=90, size=3)
ggsave("FetalSex_top50_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.detected, features = unique.top2, raster = TRUE, slot="counts", angle=90, size=3)
ggsave("FetalSex_top50_markers_counts_heatmap_rasterTRUE.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.detected, features = unique.top2, raster = TRUE, angle=90, size=3)
ggsave("FetalSex_top50_markers_logfc_heatmap_rasterTRUE.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()

Idents(seurat.object.detected) <- "Type"
cluster.averages <- AverageExpression(seurat.object.detected, assays="Spatial", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
library.averages.heatmap
ggsave("Type_unique.top2_mc_heatmap_FetalSex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("Type_unique.top2_FC_heatmap_FetalSex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
male.gene <- "DDX3Y"
top2 <- union(top2, male.gene)
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.Top50_Counts-MeanCounts_heatmap.pdf", plot = library.averages.heatmap, device = "pdf", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("FetalSex.Top50_FoldChange-MeanCounts_heatmap.pdf", plot = library.averages.heatmap, device = "pdf", width = 8, height = 10, units = "in")




setwd("..")


dir.create("DE_FetalSex_SARS-CoV-2")
setwd("DE_FetalSex_SARS-CoV-2")


Idents(seurat.object) <- "Status"
seurat.object.detected <- subset(x = seurat.object, idents = c("SARS-CoV-2"))
Idents(seurat.object.detected) <- "FetalSex"
levels(seurat.object.detected)
unfiltered.count <- table(Idents(seurat.object.detected))
unfiltered.count
write.table(unfiltered.count, "FetalSex_cellcounts.txt", sep="\t")

Idents(seurat.object.detected) <- "FetalSex"
DefaultAssay(seurat.object.detected) <- "Spatial"
seurat.object.detected <- NormalizeData(seurat.object.detected, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object.detected <- ScaleData(seurat.object.detected, features = all.genes, assay="Spatial")

dimorphism.list <- c("FCGRT", "FCGR1A", "FCGR3B","IFI6", "CXCL10", "OAS1", "CCL2", "MX1", "IL10", "IFNA1", "IFNA2" ,"IFNG", "CCL4", "CD163")
Idents(seurat.object.detected) <- "FetalSex"
cluster.averages <- AverageExpression(seurat.object.detected, assays="Spatial", slot="counts", return.seurat = TRUE, features=dimorphism.list)
write.csv(cluster.averages@assays[["Spatial"]]@counts, file = "dimorphism.list_meancounts.csv")
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = dimorphism.list, raster = FALSE, slot="counts", assay="Spatial")
ggsave("dimorphism.list_mc_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = dimorphism.list, raster = FALSE, assay="Spatial")
ggsave("dimorphism.list_FC_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.detected, features = dimorphism.list, raster = FALSE, assay="Spatial", slot="counts")
ggsave("dimorphism.list_counts_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.detected, features = dimorphism.list, raster = FALSE, assay="Spatial")
ggsave("dimorphism.list_fold-change_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

de_markers <- FindAllMarkers(seurat.object.detected, features = dimorphism.list, assay = "Spatial", only.pos = F, min.pct = 0, logfc.threshold =  0)
write.table(de_markers, "DEGs_by_FetalSex_dimorphism.list.txt", sep="\t")

de_markers <- FindAllMarkers(seurat.object.detected, features = all.genes, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_FetalSex_pos-0.693lnFC.txt", sep="\t")
## Make a mean counts object 
cluster.averages <- AverageExpression(seurat.object.detected, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
## plot the top genes
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.detected, features = unique.top2)
feature.plot
png(file=paste0("FetalSex_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.detected, features = unique.top2)
png(file=paste0("FetalSex_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
male.gene <- "DDX3Y"
top2 <- union(top2, male.gene)
unique.top2 <- unique(top2)

unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.detected, features = unique.top2)
png(file=paste0("FetalSex_top50.marker-DotPlot.png"),res=300, width=3500, height=1500)

feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="FetalSex")
dev.off()
top5.heatmap <- DoHeatmap(seurat.object.detected, features = unique.top2, raster = FALSE, angle=90, size=3)
ggsave("FetalSex_unique.top2_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.detected, features = unique.top2, raster = FALSE, slot="counts", angle=90, size=3)
ggsave("FetalSex_unique.top2_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()



Idents(seurat.object.detected) <- "Type"
cluster.averages <- AverageExpression(seurat.object.detected, assays="Spatial", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
library.averages.heatmap
ggsave("Type_unique.top2_mc_heatmap_FetalSex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("Type_unique.top2_FC_heatmap_FetalSex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
top5.heatmap <- DoHeatmap(seurat.object.detected, features = unique.top2, raster = FALSE, angle=90, size=3)
ggsave("Type_unique.top2_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()
top5.heatmap <- DoHeatmap(seurat.object.detected, features = unique.top2, raster = FALSE, slot="counts", angle=90, size=3)
ggsave("Type_unique.top2_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 8, height = 10, units = "in") + NoLegend()


Idents(seurat.object.detected) <- "FetalSex"



val.list <- c("TIMP1", "PRG2", "AOC1", "NOTUM","RPS4Y1")


val.list <- c("TIMP1", "PRG2", "AOC1", "NOTUM","RPS4Y1")

## 	REDO SPATIAL PLOTS ABOVE
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
image.list <- c("S01")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S01.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S03.1")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S03.1.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S04.2")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S04.2.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S15.3")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S15.3.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")

image.list <- c("S16.4")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S16.4.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S17.5")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S17.5.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S18.6")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S18.6.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S19.7")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S19.7.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S20.8")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S20.8.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S21.9")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S21.9.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")
image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")

image.list <- c("S21.9")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S21.9.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S23a.11")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S23a.11.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S23b.12")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S23b.12.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S24.13")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S24.13.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S25.14")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S25.14.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")

image.list <- c("S26.15")
p1<- SpatialPlot(seurat.object, features = val.list, images = image.list) + patchwork::plot_layout(ncol = 2)
p1
ggsave("val.list_S26.15.pdf", plot = p1, device = "pdf", width = 8, height = 11, units = "in")





# seurat.object.detected <- ScaleData(seurat.object.detected, features = all.genes, assay="Spatial")
Idents(seurat.object.detected) <- "FetalSex"
cluster.averages <- AverageExpression(seurat.object.detected, assays="Spatial", slot="counts", return.seurat = TRUE, features=val.list)
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = val.list, raster = FALSE, slot="counts", assay="Spatial")
library.averages.heatmap
ggsave("val.list_mc_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = val.list, raster = FALSE, assay="Spatial")
ggsave("val.list_FC_heatmap_fetalsex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.detected, features = val.list, raster = FALSE, assay="Spatial", slot="counts")
ggsave("val.list_counts_heatmap_fetalsex.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.detected, features = val.list, raster = FALSE, assay="Spatial")
ggsave("val.list_fold-change_heatmap_fetalsex.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")

Idents(seurat.object.detected) <- "Type"
cluster.averages <- AverageExpression(seurat.object.detected, assays="Spatial", slot="counts", return.seurat = TRUE, features=val.list)
####Use this to determine the sum of mean counts for viral transcripts per cluster and rank viral transcripts based on overall viral transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = val.list, raster = FALSE, slot="counts", assay="Spatial")
library.averages.heatmap
ggsave("val.list_mc_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for viral transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = val.list, raster = FALSE, assay="Spatial")
ggsave("val.list_FC_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.detected, features = val.list, raster = FALSE, assay="Spatial", slot="counts")
ggsave("val.list_counts_heatmap_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object.detected, features = val.list, raster = FALSE, assay="Spatial")
ggsave("val.list_fold-change_heatmap_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 11, height = 8, units = "in")




image.list <- c("S01", "S03.1", "S04.2", "S15.3", "S16.4", "S17.5", "S18.6", "S19.7", "S20.8", "S21.9", "S22.10", "S23a.11", "S23b.12", "S24.13", "S25.14", "S26.15")
p5 <- SpatialFeaturePlot(seurat.object, features = "TIMP1", images=image.list, ncol=3, min.cutoff=0, max.cutoff=7) + patchwork::plot_layout(ncol = 4)
ggsave("TIMP1_SpatialFeatureplot-cropped-scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "PRG2", images=image.list, ncol=3, min.cutoff=0, max.cutoff=7) + patchwork::plot_layout(ncol = 4)
ggsave("PRG2_SpatialFeatureplot-cropped-scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "AOC1", images=image.list, ncol=3, min.cutoff=0, max.cutoff=7) + patchwork::plot_layout(ncol = 4)
ggsave("AOC1_SpatialFeatureplot-cropped-scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "NOTUM", images=image.list, ncol=3, min.cutoff=0, max.cutoff=7) + patchwork::plot_layout(ncol = 4)
ggsave("NOTUM_SpatialFeatureplot-cropped-scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")
p5 <- SpatialFeaturePlot(seurat.object, features = "RPS4Y1", images=image.list, ncol=3, min.cutoff=0, max.cutoff=7) + patchwork::plot_layout(ncol = 4)
ggsave("RPS4Y1_SpatialFeatureplot-cropped-scaleyfixey.pdf", plot = p5, device = "pdf", width = 9, height = 10, units = "in")





## Make a means count heatmap using cluster.averages
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.Top10_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="Spatial")
ggsave("FetalSex.Top10.2_FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="Spatial")
ggsave("FetalSex.SARS-CoV-2_Counts-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="Spatial")
ggsave("FetalSex.SARS-CoV-2._FoldChange-MeanCounts_heatmap.png", plot = library.averages.heatmap, device = "png", width = 8, height = 10, units = "in")

de_markers <- FindAllMarkers(seurat.object.detected, features = all.genes, assay = "Spatial", only.pos = F, min.pct = 0.10, logfc.threshold =  0.2)
write.table(de_markers, "DEGs_by_FetalSex_0.2lnFC.txt", sep="\t")







setwd("..")




######################## ######################## ########################
######################## ######################## ########################
######################## ######################## ########################

set.seed(seed=1)
setwd("/home/ebarrozo/NCS/results")
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
library(SeuratDisk)

######################## ######################## ########################
######################## Had issues installing monocle3 on aagaardlab3, so saved csf-macrophage-subset and doing monocle3 on local
setwd("/Users/enricobarrozo/Library/CloudStorage/Box-Box/AagaardLab/Visium/SpaceRanger/results/seurat_human_v2/subset.macs/")
	## start with cluster 9 based on HSC markers


library(ggplot2)
library(Seurat)
load("human_spatial_data-macs-annotated_v1.RData")


## run immunology on polarization groups and clusters
dir.create("immunology.clusters.Polarization")
setwd("immunology.clusters.Polarization")

dir.create("immunology.clusters.clusters")
setwd("immunology.clusters.clusters")


Idents(seurat.object) <- "seurat_clusters"
levels(seurat.object)
#  [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12"


Idents(seurat.object) <- "Polarization"
levels(seurat.object)
	# [1] "M2"    "M0.M1" "M2a"   "M0.M2" "M1"    "M0.1"  "M0"    "DC"    "M2.M1"

## use RNA for all DE analysis and plots
# seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "SCT")
seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
DefaultAssay(seurat.object) <- "Spatial"
library("Nebulosa")


## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("ITGAM", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_ITGAM_THBD.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object, "ITGAM")
p2 <- FeaturePlot(seurat.object, features="ITGAM", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object, "THBD")
p2 <- FeaturePlot(seurat.object, features="THBD", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_THBD.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD4", "IFNG"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD4_IFNG.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object, "CD4")
p2 <- FeaturePlot(seurat.object, features="CD4", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD4.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object, "IFNG")
p2 <- FeaturePlot(seurat.object, features="IFNG", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_IFNG.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD38", "PTPRC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD38_PTPRC.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object, "CD38")
p2 <- FeaturePlot(seurat.object, features="CD38", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD38.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object, "PTPRC")
p2 <- FeaturePlot(seurat.object, features="PTPRC", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_PTPRC.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD38_FPR2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object, "CD38")
p2 <- FeaturePlot(seurat.object, features="CD38", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD38.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object, "FPR2")
p2 <- FeaturePlot(seurat.object, features="FPR2", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_FPR2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")



## kernel density estimation 
# https://www.nature.com/articles/s41591-021-01329-2#Sec8
major.t <- c("CD4", "CD8A", "CCR7", "PTPRC", "SELL", "CD27", "CD38", "CD44", "CXCR5", "CD40LG", "CCR7", "FOXP3", "IKZF2")
gd.t<- c("TRGV9", "TRDV2")
MAIT.t <- c("TRAV10", "TRAJ18", "TRAV1-2")
NK.t <- c("CD3G", "FCGR3A", "NCAM1", "NCR1")
TH1.t <- c("IFNG", "TBX21", "TNF")
TH2.t <- c("GATA3", "IL4")
TH17.t <- c("RORC", "IL17A", "IL17F", "IL21")
## https://www.bioconductor.org/packages/release/bioc/vignettes/Nebulosa/inst/doc/nebulosa_seurat.html#1_Overview
  ## https://academic.oup.com/bioinformatics/article-abstract/37/16/2485/6103785
# BiocManager::install("Nebulosa")
library("Nebulosa")

# SARS.genes <- c("ORF1AB-1-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
p1 <- plot_density(seurat.object, "CD4")
p2 <- FeaturePlot(seurat.object, "CD8A")
p5 <- CombinePlots(plots = list(p1, p2))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("UMAP_KernelDensity_CD4-CD8A.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots
p3 <- plot_density(seurat.object, c("CD4", "CD8A"))
p5 <- p3 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_CD4-CD8A_Joint.pdf", plot = p5, device = "pdf", width = 5, height = 7, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD4", "CD8A"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_major.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("TRGV9", "TRDV2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_gd.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("TRAV10", "TRAJ18"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MAIT.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD3G", "FCGR3A","NCAM1", "NCR1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_NK.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")

	#	Error in linbin2D.ks(x, gpoints[[1]], gpoints[[2]], w = w) : 
	#	  NA/NaN/Inf in foreign function call (arg 10)
	#	In addition: Warning message:
	#	In ks.defaults(x = x, w = w, binned = binned, bgridsize = bgridsize,  :
	#	  Weights don't sum to sample size - they have been scaled accordingly


## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("IFNG", "TBX21", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH1.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 16, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("GATA3", "IL4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH2.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("RORC", "IL17A", "IL17F", "IL21"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH17.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")




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
# seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

############## Activation gating strategy https://www.nature.com/articles/s41590-021-01049-2/figures/10
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD3D", "NCAM1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TorNK_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD4", "KLRB1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRB1-CD161.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD4", "CD38"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD38.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object, c("CD4", "CD69"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD69.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD4", "KLRK1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRK1-NKG2D.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD4", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("FOXP3", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-FOXP3-CD4.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("KIT", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("TNF", "TNFRSF8"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.TNFRSF8.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")


p4 <- plot_density(seurat.object, c("CD14", "ITGAM"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("SIGLEC1", "FOLR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("SIGLEC1", "DDX3Y"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC-fetalsex_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
	## Placenta-associated maternal macrophages PAMMS (FOLR2 neg, HLA positive); HBCs FOLR2+HLA- ; Thomas et al. Naomi McGovern, JEM 2020
p4 <- plot_density(seurat.object, c("FOLR2", "HLA-DRA"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMMs_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM1_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM2_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("HLA-DRA", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM1-HLA_CCR2_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")



p4 <- plot_density(seurat.object, c("ITGAM", "CD14"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M0.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD80", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD11B", "MYC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD163", "MRC1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2a.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("HLA-DRA", "IL6"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD14", "CD163"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2c.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD14", "CD1C"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object, c("CD14", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p1 <- plot_density(seurat.object, "THBD")
p2 <- FeaturePlot(seurat.object, "THBD")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "ITGB3")
p2 <- FeaturePlot(seurat.object, "ITGB3")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_Megakaryocyte.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "CD1C")
p2 <- FeaturePlot(seurat.object, "CD1C")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "IL6")
p2 <- FeaturePlot(seurat.object, "IL6")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "TNF")
p2 <- FeaturePlot(seurat.object, "TNF")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "MYC")
p2 <- FeaturePlot(seurat.object, "MYC")
p5 <- CombinePlots(plots = list(p1, p2))


ggsave("UMAP_M2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "FOLR2")
p2 <- FeaturePlot(seurat.object, "FOLR2")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_HBCs-FOLR2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "HLA-DRA")
p2 <- FeaturePlot(seurat.object, "HLA-DRA")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_PAMM-HLA-DRA.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "CD33")
p2 <- FeaturePlot(seurat.object, "CD33")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeliodProgenitor-M0.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p4 <- plot_density(seurat.object, c("KIT", "CD33"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_MyeliodProgenitor-M0-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object, "CCR2")
p2 <- FeaturePlot(seurat.object, "CCR2")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_PAMM-CCR2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object, "ITGAM")
p2 <- FeaturePlot(seurat.object, "ITGAM")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_PAMM-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "CD14")
p2 <- FeaturePlot(seurat.object, "CD14")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_PAMM-CD14.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object, "RNASE1")
p2 <- FeaturePlot(seurat.object, "RNASE1")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_YS_Mac-RNASE1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "HBE1")
p2 <- FeaturePlot(seurat.object, "HBE1")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_YS_Mac1-HBE1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "PLOD2")
p2 <- FeaturePlot(seurat.object, "PLOD2")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_YS_Mac1-PLOD2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "BNIP3")
p2 <- FeaturePlot(seurat.object, "BNIP3")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_YS_Mac1-BNIP3.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "GSTA1")
p2 <- FeaturePlot(seurat.object, "GSTA1")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_YS_Mac2-GSTA1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "FGA")
p2 <- FeaturePlot(seurat.object, "FGA")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_YS_Mac2-FGA.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object, "FGG")
p2 <- FeaturePlot(seurat.object, "FGG")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_YS_Mac2-FGG.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")



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


M2.list<- c("MYC", "EGR2", "CD163", "MRC1", "TGFB1", "VEGFA", "HLA-DRA", "IL6")
seurat.object <- AddModuleScore(object = seurat.object, features = M2.list, ctrl = 5, name = 'M2.list.score')
p1 <- FeaturePlot(seurat.object, features = 'M2.list.score1')
ggsave("FeaturePlot_M2.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'M2.list.score1')
ggsave("RidgePlot_M2.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "M2.list.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M2.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "M2.list.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M2.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")

M1.list<- c("CD38", "FPR2")
seurat.object <- AddModuleScore(object = seurat.object, features = M1.list, ctrl = 5, name = 'M1.list.score')
p1 <- FeaturePlot(seurat.object, features = 'M1.list.score1')
ggsave("FeaturePlot_M1.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1<-RidgePlot(seurat.object, features = 'M1.list.score1')
ggsave("RidgePlot_M1.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- VlnPlot(seurat.object, features = "M1.list.score1", pt.size = 0.1) + NoLegend()
ggsave("VlnPlot_M1.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")
p1 <- FeatureScatter(seurat.object, feature1 = "M1.list.score1", feature2 = "percent.viral")
ggsave("FeatureScatter_M1.list.png", plot = p1, device = "png", width = 7, height = 5, units = "in")


setwd("..")
