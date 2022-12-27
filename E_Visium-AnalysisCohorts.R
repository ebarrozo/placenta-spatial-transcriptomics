## covid-visium_revisions_v2.R

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
## revisions_v2: Rename SampleIDs to SampleCodes, rename condition to AnalysisCohort
## Analyze data on AagaardLab3
# 


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
dir.create("revisions_v2")
setwd("revisions_v2")

library(ggpubr)
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

setwd("/home/ebarrozo/visium/results/seurat_human_v2")
load("/home/ebarrozo/visium/results/seurat_human_v2/annotated/human_spatial_data-integrated-annotated_v1.RData")
DefaultAssay(seurat.object) <- "Spatial"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object) #  "S01"  "S03"  "S04"  "S15"  "S16"  "S17"  "S18"  "S19"  "S20"  "S21"  "S22" "S23a" "S23b" "S24"  "S25"  "S26"
new.metadata <- c("F", "F", "F", "F", 
	"M", "F", "M", "M", 
	"F", "F", "M", "F", 
	"F", "M","M","M")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$FetalSex <- Idents(seurat.object)
Idents(seurat.object) <- "FetalSex"


setwd("/home/ebarrozo/visium/results/seurat_human_v2/revisions_v2")

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("NC1a", "NC1b", "NC1c", "ND1", "SP4", "SP1", "SP2", "SP3", "ND2", "ND3", "SP5", "HP1a", "HP1b", "NC2","NC3","NC4")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$SampleCode <- Idents(seurat.object)
Idents(seurat.object) <- "SampleCode"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("NegativeControl", "NegativeControl", "NegativeControl", "NotDetected", "SparsePositive", "SparsePositive", "SparsePositive", "SparsePositive", "NotDetected", "NotDetected", "SparsePositive", "HighPositive", "HighPositive", "NegativeControl","NegativeControl","NegativeControl")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$AnalysisCohort <- Idents(seurat.object)
Idents(seurat.object) <- "AnalysisCohort"

p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "SampleCode",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='SampleCode')
p2
ggsave("UMAP-SampleCode.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
ggsave("UMAP-SampleCodeXL.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")

p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "AnalysisCohort",  raster=FALSE, cols = mypal3) + labs(title = NULL, color='AnalysisCohort')
p2
ggsave("UMAP-AnalysisCohort.pdf", plot = p2, device = "pdf", width = 8.5, height = 6, units = "in")
ggsave("UMAP-AnalysisCohort_XL.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")

p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Subtype", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Subtype')
ggsave("UMAP-Subtype.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")


dir.create("immunology.clusters_clusters")
setwd("immunology.clusters_clusters")
Idents(seurat.object) <- "seurat_clusters"
## use RNA for all DE analysis and plots
# seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "SCT")
seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
DefaultAssay(seurat.object) <- "Spatial"


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
Idents(seurat.object) <- "seurat_clusters"

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

setwd("..")


dir.create("immunology.clusters_Type")
setwd("immunology.clusters_Type")
Idents(seurat.object) <- "Type"
## use RNA for all DE analysis and plots
# seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "SCT")
# seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
DefaultAssay(seurat.object) <- "Spatial"


Idents(seurat.object) <- "Type"
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

Idents(seurat.object) <- "Type"
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

Idents(seurat.object) <- "Type"
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

Idents(seurat.object) <- "Type"
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
Idents(seurat.object) <- "Type"

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

setwd("..")

############################################################# ############################################################# #############################################################
############################################################# DE Analysis Btwn AnalysisCohort #############################################################
############################################################# ############################################################# #############################################################

dir.create("DE_AnalysisCohort")
setwd("DE_AnalysisCohort")
Idents(seurat.object)
Idents(seurat.object) <- "AnalysisCohort"

p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "AnalysisCohort")
ggsave("integrated_UMAP-AnalysisCohort.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")
p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Subtype", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Subtype')
ggsave("UMAP-Subtype.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p4 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "SampleCode", cols = mypal3, raster=FALSE) + labs(title = NULL, color='SampleCode')
ggsave("UMAP-SampleCode.pdf".subset, plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = all.genes, assay="Spatial")

Idents(seurat.object) <- "AnalysisCohort"
de_markers <- FindAllMarkers(seurat.object, features = all.genes, assay = "Spatial", only.pos = FALSE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "DEGs_byAnalysisCohort_2log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = FALSE)
ggsave("top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = FALSE)
ggsave("top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("top10.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="AnalysisCohort")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("top15.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="AnalysisCohort")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="AnalysisCohort")
dev.off()



setwd("..")

############################################################# ############################################################# #############################################################
############################################################# DE Analysis Btwn AnalysisCohort #############################################################
############################################################# ############################################################# #############################################################

dir.create("DE_AnalysisCohort-SparsePositiveVNegativeControl")
setwd("DE_AnalysisCohort-SparsePositiveVNegativeControl")
Idents(seurat.object)
Idents(seurat.object) <- "AnalysisCohort"

seurat.object.subset <- subset(x = seurat.object, idents = c("SparsePositive", "NegativeControl"))


p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "AnalysisCohort")
ggsave("integrated_UMAP-AnalysisCohort.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")
p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Subtype", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Subtype')
ggsave("UMAP-Subtype.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p4 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "SampleCode", cols = mypal3, raster=FALSE) + labs(title = NULL, color='SampleCode')
ggsave("UMAP-SampleCode.pdf".subset, plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object.subset) <- "Spatial"
seurat.object.subset <- NormalizeData(seurat.object.subset, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object.subset <- ScaleData(seurat.object.subset, features = all.genes, assay="Spatial")

Idents(seurat.object.subset) <- "AnalysisCohort"
de_markers <- FindAllMarkers(seurat.object.subset, features = all.genes, assay = "Spatial", only.pos = FALSE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "DEGs_byAnalysisCohort_2log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, raster = FALSE)
ggsave("top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, slot="counts", raster = FALSE)
ggsave("top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top10.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="AnalysisCohort")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top15.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="AnalysisCohort")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="AnalysisCohort")
dev.off()



setwd("..")
############################################################# ############################################################# #############################################################
############################################################# DE Analysis Btwn AnalysisCohort #############################################################
############################################################# ############################################################# #############################################################

dir.create("DE_AnalysisCohort-NotDetectedVNegativeControl")
setwd("DE_AnalysisCohort-NotDetectedVNegativeControl")
Idents(seurat.object)
Idents(seurat.object) <- "AnalysisCohort"

seurat.object.subset <- subset(x = seurat.object, idents = c("NotDetected", "NegativeControl"))


p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "AnalysisCohort")
ggsave("integrated_UMAP-AnalysisCohort.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")
p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Subtype", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Subtype')
ggsave("UMAP-Subtype.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p4 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "SampleCode", cols = mypal3, raster=FALSE) + labs(title = NULL, color='SampleCode')
ggsave("UMAP-SampleCode.pdf".subset, plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object.subset) <- "Spatial"
seurat.object.subset <- NormalizeData(seurat.object.subset, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object.subset <- ScaleData(seurat.object.subset, features = all.genes, assay="Spatial")

Idents(seurat.object.subset) <- "AnalysisCohort"
de_markers <- FindAllMarkers(seurat.object.subset, features = all.genes, assay = "Spatial", only.pos = FALSE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "DEGs_byAnalysisCohort_2log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, raster = FALSE)
ggsave("top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, slot="counts", raster = FALSE)
ggsave("top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top10.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="AnalysisCohort")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top15.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="AnalysisCohort")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="AnalysisCohort")
dev.off()



setwd("..")

############################################################# ############################################################# #############################################################
############################################################# DE Analysis Btwn AnalysisCohort #############################################################
############################################################# ############################################################# #############################################################

dir.create("DE_AnalysisCohort-HighPositiveVNegativeControl")
setwd("DE_AnalysisCohort-HighPositiveVNegativeControl")
Idents(seurat.object)
Idents(seurat.object) <- "AnalysisCohort"

seurat.object.subset <- subset(x = seurat.object, idents = c("HighPositive", "NegativeControl"))


p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Subtype", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Subtype')
ggsave("UMAP-Subtype.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p4 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "SampleCode", cols = mypal3, raster=FALSE) + labs(title = NULL, color='SampleCode')
ggsave("UMAP-SampleCode.pdf".subset, plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")


p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "AnalysisCohort")
ggsave("UMAP-AnalysisCohort.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object.subset) <- "Spatial"
seurat.object.subset <- NormalizeData(seurat.object.subset, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object.subset <- ScaleData(seurat.object.subset, features = all.genes, assay="Spatial")

Idents(seurat.object.subset) <- "AnalysisCohort"
de_markers <- FindAllMarkers(seurat.object.subset, features = all.genes, assay = "Spatial", only.pos = FALSE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "DEGs_byAnalysisCohort_2log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, raster = FALSE)
ggsave("top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, slot="counts", raster = FALSE)
ggsave("top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top10.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="AnalysisCohort")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top15.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="AnalysisCohort")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="AnalysisCohort")
dev.off()



setwd("..")



############################################################# ############################################################# #############################################################
############################################################# Subset analysis of HighPositive #############################################################
############################################################# ############################################################# #############################################################


dir.create("Subset-HighPositive")
setwd("Subset-HighPositive")
Idents(seurat.object)
Idents(seurat.object) <- "AnalysisCohort"

seurat.object.subset <- subset(x = seurat.object, idents = c("HighPositive"))


p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Subtype", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Subtype')
ggsave("UMAP-Subtype.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p4 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "SampleCode", cols = mypal3, raster=FALSE) + labs(title = NULL, color='SampleCode')
ggsave("UMAP-SampleCode.pdf".subset, plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")
p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "AnalysisCohort")
ggsave("UMAP-AnalysisCohort.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

Idents(seurat.object.subset) <- "Type"
de_markers <- FindAllMarkers(seurat.object.subset, features = all.genes, assay = "Spatial", only.pos = FALSE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "DEGs_byType_2log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, raster = FALSE)
ggsave("top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, slot="counts", raster = FALSE)
ggsave("top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top10.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top15.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type")
dev.off()
Idents(seurat.object.subset) <- "Type"
## use RNA for all DE analysis and plots
# seurat.object.subset <- ScaleData(seurat.object.subset, features = all.genes, assay = "SCT")
# seurat.object.subset <- ScaleData(seurat.object.subset, features = all.genes, assay = "Spatial")
DefaultAssay(seurat.object.subset) <- "Spatial"
Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("ITGAM", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_ITGAM_THBD.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "ITGAM")
p2 <- FeaturePlot(seurat.object.subset, features="ITGAM", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "THBD")
p2 <- FeaturePlot(seurat.object.subset, features="THBD", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_THBD.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "IFNG"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD4_IFNG.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "CD4")
p2 <- FeaturePlot(seurat.object.subset, features="CD4", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD4.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "IFNG")
p2 <- FeaturePlot(seurat.object.subset, features="IFNG", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_IFNG.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD38", "PTPRC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD38_PTPRC.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "CD38")
p2 <- FeaturePlot(seurat.object.subset, features="CD38", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD38.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "PTPRC")
p2 <- FeaturePlot(seurat.object.subset, features="PTPRC", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_PTPRC.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD38_FPR2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "CD38")
p2 <- FeaturePlot(seurat.object.subset, features="CD38", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD38.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "FPR2")
p2 <- FeaturePlot(seurat.object.subset, features="FPR2", slot="scale.data")
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
Idents(seurat.object.subset) <- "Type"

# SARS.genes <- c("ORF1AB-1-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
p1 <- plot_density(seurat.object.subset, "CD4")
p2 <- FeaturePlot(seurat.object.subset, "CD8A")
p5 <- CombinePlots(plots = list(p1, p2))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("UMAP_KernelDensity_CD4-CD8A.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots
p3 <- plot_density(seurat.object.subset, c("CD4", "CD8A"))
p5 <- p3 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_CD4-CD8A_Joint.pdf", plot = p5, device = "pdf", width = 5, height = 7, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "CD8A"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_major.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("TRGV9", "TRDV2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_gd.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("TRAV10", "TRAJ18"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MAIT.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD3G", "FCGR3A","NCAM1", "NCR1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_NK.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")

	#	Error in linbin2D.ks(x, gpoints[[1]], gpoints[[2]], w = w) : 
	#	  NA/NaN/Inf in foreign function call (arg 10)
	#	In addition: Warning message:
	#	In ks.defaults(x = x, w = w, binned = binned, bgridsize = bgridsize,  :
	#	  Weights don't sum to sample size - they have been scaled accordingly


## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("IFNG", "TBX21", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH1.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 16, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("GATA3", "IL4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH2.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("RORC", "IL17A", "IL17F", "IL21"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH17.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")




## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(seurat.object.subset, features = lymphoid.lineage)
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
feature.plot <- DotPlot(seurat.object.subset, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(seurat.object.subset, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

macrophage.lineage <- c("CD14", "CD68", "CSF1R", "ADGRE1", "SIGLEC1", "CD38", "GPR18", "FPR2", "MYC", "EGR2")
feature.plot <- DotPlot(seurat.object.subset, features = macrophage.lineage)
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

feature.plot <- DotPlot(seurat.object.subset, features = macrophage.lineage.markers.2)
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

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=t.lineages)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, slot="counts", assay="Spatial")
ggsave("t.lineages_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, assay="Spatial")
ggsave("t.lineages_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=macrophage.lineage.markers.2)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("macrophage.lineage.markers.2_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, assay="Spatial")
ggsave("macrophage.lineage.markers.2_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=lymphoid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("lymphoid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, assay="Spatial")
ggsave("lymphoid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=myeloid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("myeloid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, assay="Spatial")
ggsave("myeloid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
# seurat.object.subset <- DietSeurat(seurat.object.subset, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

############## Activation gating strategy https://www.nature.com/articles/s41590-021-01049-2/figures/10
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD3D", "NCAM1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TorNK_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "KLRB1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRB1-CD161.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "CD38"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD38.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "CD69"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD69.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD4", "KLRK1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRK1-NKG2D.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD4", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("FOXP3", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-FOXP3-CD4.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("KIT", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("TNF", "TNFRSF8"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.TNFRSF8.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")


p4 <- plot_density(seurat.object.subset, c("CD14", "ITGAM"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("SIGLEC1", "FOLR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("SIGLEC1", "DDX3Y"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC-fetalsex_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
	## Placenta-associated maternal macrophages PAMMS (FOLR2 neg, HLA positive); HBCs FOLR2+HLA- ; Thomas et al. Naomi McGovern, JEM 2020
p4 <- plot_density(seurat.object.subset, c("FOLR2", "HLA-DRA"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMMs_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM1_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM2_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("ITGAM", "CD14"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M0.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD80", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD11B", "MYC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD163", "MRC1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2a.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("HLA-DRA", "IL6"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD14", "CD163"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2c.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD14", "CD1C"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD14", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p1 <- plot_density(seurat.object.subset, "THBD")
p2 <- FeaturePlot(seurat.object.subset, "THBD")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "ITGB3")
p2 <- FeaturePlot(seurat.object.subset, "ITGB3")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_Megakaryocyte.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "CD1C")
p2 <- FeaturePlot(seurat.object.subset, "CD1C")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "IL6")
p2 <- FeaturePlot(seurat.object.subset, "IL6")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "TNF")
p2 <- FeaturePlot(seurat.object.subset, "TNF")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "MYC")
p2 <- FeaturePlot(seurat.object.subset, "MYC")
p5 <- CombinePlots(plots = list(p1, p2))


ggsave("UMAP_M2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "FOLR2")
p2 <- FeaturePlot(seurat.object.subset, "FOLR2")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_HBCs-FOLR2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "HLA-DRA")
p2 <- FeaturePlot(seurat.object.subset, "HLA-DRA")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_PAMM-HLA-DRA.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "CD33")
p2 <- FeaturePlot(seurat.object.subset, "CD33")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeliodProgenitor-M0.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p4 <- plot_density(seurat.object.subset, c("KIT", "CD33"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_MyeliodProgenitor-M0-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

setwd("..")


############################################################# ############################################################# #############################################################
############################################################# Subset analysis of NotDetected #############################################################
############################################################# ############################################################# #############################################################


dir.create("Subset-NotDetected")
setwd("Subset-NotDetected")
Idents(seurat.object)
Idents(seurat.object) <- "AnalysisCohort"

seurat.object.subset <- subset(x = seurat.object, idents = c("NotDetected"))


p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Subtype", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Subtype')
ggsave("UMAP-Subtype.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p4 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "SampleCode", cols = mypal3, raster=FALSE) + labs(title = NULL, color='SampleCode')
ggsave("UMAP-SampleCode.pdf".subset, plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")
p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "AnalysisCohort")
ggsave("UMAP-AnalysisCohort.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

Idents(seurat.object.subset) <- "Type"
de_markers <- FindAllMarkers(seurat.object.subset, features = all.genes, assay = "Spatial", only.pos = FALSE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "DEGs_byType_2log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, raster = FALSE)
ggsave("top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, slot="counts", raster = FALSE)
ggsave("top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top10.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top15.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type")
dev.off()
Idents(seurat.object.subset) <- "Type"
## use RNA for all DE analysis and plots
# seurat.object.subset <- ScaleData(seurat.object.subset, features = all.genes, assay = "SCT")
# seurat.object.subset <- ScaleData(seurat.object.subset, features = all.genes, assay = "Spatial")
DefaultAssay(seurat.object.subset) <- "Spatial"
Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("ITGAM", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_ITGAM_THBD.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "ITGAM")
p2 <- FeaturePlot(seurat.object.subset, features="ITGAM", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "THBD")
p2 <- FeaturePlot(seurat.object.subset, features="THBD", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_THBD.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "IFNG"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD4_IFNG.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "CD4")
p2 <- FeaturePlot(seurat.object.subset, features="CD4", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD4.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "IFNG")
p2 <- FeaturePlot(seurat.object.subset, features="IFNG", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_IFNG.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD38", "PTPRC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD38_PTPRC.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "CD38")
p2 <- FeaturePlot(seurat.object.subset, features="CD38", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD38.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "PTPRC")
p2 <- FeaturePlot(seurat.object.subset, features="PTPRC", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_PTPRC.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD38_FPR2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "CD38")
p2 <- FeaturePlot(seurat.object.subset, features="CD38", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD38.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "FPR2")
p2 <- FeaturePlot(seurat.object.subset, features="FPR2", slot="scale.data")
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
Idents(seurat.object.subset) <- "Type"

# SARS.genes <- c("ORF1AB-1-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
p1 <- plot_density(seurat.object.subset, "CD4")
p2 <- FeaturePlot(seurat.object.subset, "CD8A")
p5 <- CombinePlots(plots = list(p1, p2))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("UMAP_KernelDensity_CD4-CD8A.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots
p3 <- plot_density(seurat.object.subset, c("CD4", "CD8A"))
p5 <- p3 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_CD4-CD8A_Joint.pdf", plot = p5, device = "pdf", width = 5, height = 7, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "CD8A"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_major.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("TRGV9", "TRDV2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_gd.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("TRAV10", "TRAJ18"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MAIT.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD3G", "FCGR3A","NCAM1", "NCR1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_NK.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")

	#	Error in linbin2D.ks(x, gpoints[[1]], gpoints[[2]], w = w) : 
	#	  NA/NaN/Inf in foreign function call (arg 10)
	#	In addition: Warning message:
	#	In ks.defaults(x = x, w = w, binned = binned, bgridsize = bgridsize,  :
	#	  Weights don't sum to sample size - they have been scaled accordingly


## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("IFNG", "TBX21", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH1.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 16, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("GATA3", "IL4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH2.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("RORC", "IL17A", "IL17F", "IL21"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH17.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")




## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(seurat.object.subset, features = lymphoid.lineage)
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
feature.plot <- DotPlot(seurat.object.subset, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(seurat.object.subset, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

macrophage.lineage <- c("CD14", "CD68", "CSF1R", "ADGRE1", "SIGLEC1", "CD38", "GPR18", "FPR2", "MYC", "EGR2")
feature.plot <- DotPlot(seurat.object.subset, features = macrophage.lineage)
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

feature.plot <- DotPlot(seurat.object.subset, features = macrophage.lineage.markers.2)
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

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=t.lineages)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, slot="counts", assay="Spatial")
ggsave("t.lineages_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, assay="Spatial")
ggsave("t.lineages_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=macrophage.lineage.markers.2)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("macrophage.lineage.markers.2_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, assay="Spatial")
ggsave("macrophage.lineage.markers.2_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=lymphoid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("lymphoid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, assay="Spatial")
ggsave("lymphoid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=myeloid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("myeloid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, assay="Spatial")
ggsave("myeloid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
# seurat.object.subset <- DietSeurat(seurat.object.subset, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

############## Activation gating strategy https://www.nature.com/articles/s41590-021-01049-2/figures/10
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD3D", "NCAM1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TorNK_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "KLRB1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRB1-CD161.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "CD38"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD38.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "CD69"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD69.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD4", "KLRK1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRK1-NKG2D.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD4", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("FOXP3", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-FOXP3-CD4.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("KIT", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("TNF", "TNFRSF8"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.TNFRSF8.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")


p4 <- plot_density(seurat.object.subset, c("CD14", "ITGAM"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("SIGLEC1", "FOLR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("SIGLEC1", "DDX3Y"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC-fetalsex_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
	## Placenta-associated maternal macrophages PAMMS (FOLR2 neg, HLA positive); HBCs FOLR2+HLA- ; Thomas et al. Naomi McGovern, JEM 2020
p4 <- plot_density(seurat.object.subset, c("FOLR2", "HLA-DRA"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMMs_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM1_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM2_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("ITGAM", "CD14"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M0.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD80", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD11B", "MYC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD163", "MRC1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2a.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("HLA-DRA", "IL6"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD14", "CD163"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2c.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD14", "CD1C"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD14", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p1 <- plot_density(seurat.object.subset, "THBD")
p2 <- FeaturePlot(seurat.object.subset, "THBD")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "ITGB3")
p2 <- FeaturePlot(seurat.object.subset, "ITGB3")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_Megakaryocyte.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "CD1C")
p2 <- FeaturePlot(seurat.object.subset, "CD1C")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "IL6")
p2 <- FeaturePlot(seurat.object.subset, "IL6")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "TNF")
p2 <- FeaturePlot(seurat.object.subset, "TNF")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "MYC")
p2 <- FeaturePlot(seurat.object.subset, "MYC")
p5 <- CombinePlots(plots = list(p1, p2))


ggsave("UMAP_M2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "FOLR2")
p2 <- FeaturePlot(seurat.object.subset, "FOLR2")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_HBCs-FOLR2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "HLA-DRA")
p2 <- FeaturePlot(seurat.object.subset, "HLA-DRA")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_PAMM-HLA-DRA.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "CD33")
p2 <- FeaturePlot(seurat.object.subset, "CD33")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeliodProgenitor-M0.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p4 <- plot_density(seurat.object.subset, c("KIT", "CD33"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_MyeliodProgenitor-M0-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

setwd("..")


############################################################# ############################################################# #############################################################
############################################################# Subset analysis of SparsePositive #############################################################
############################################################# ############################################################# #############################################################


dir.create("Subset-SparsePositive")
setwd("Subset-SparsePositive")
Idents(seurat.object)
Idents(seurat.object) <- "AnalysisCohort"

seurat.object.subset <- subset(x = seurat.object, idents = c("SparsePositive"))


p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Subtype", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Subtype')
ggsave("UMAP-Subtype.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p4 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "SampleCode", cols = mypal3, raster=FALSE) + labs(title = NULL, color='SampleCode')
ggsave("UMAP-SampleCode.pdf".subset, plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")
p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "AnalysisCohort")
ggsave("UMAP-AnalysisCohort.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

Idents(seurat.object.subset) <- "Type"
de_markers <- FindAllMarkers(seurat.object.subset, features = all.genes, assay = "Spatial", only.pos = FALSE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "DEGs_byType_2log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, raster = FALSE)
ggsave("top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, slot="counts", raster = FALSE)
ggsave("top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top10.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top15.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type")
dev.off()
Idents(seurat.object.subset) <- "Type"
## use RNA for all DE analysis and plots
# seurat.object.subset <- ScaleData(seurat.object.subset, features = all.genes, assay = "SCT")
# seurat.object.subset <- ScaleData(seurat.object.subset, features = all.genes, assay = "Spatial")
DefaultAssay(seurat.object.subset) <- "Spatial"
Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("ITGAM", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_ITGAM_THBD.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "ITGAM")
p2 <- FeaturePlot(seurat.object.subset, features="ITGAM", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "THBD")
p2 <- FeaturePlot(seurat.object.subset, features="THBD", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_THBD.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "IFNG"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD4_IFNG.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "CD4")
p2 <- FeaturePlot(seurat.object.subset, features="CD4", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD4.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "IFNG")
p2 <- FeaturePlot(seurat.object.subset, features="IFNG", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_IFNG.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD38", "PTPRC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD38_PTPRC.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "CD38")
p2 <- FeaturePlot(seurat.object.subset, features="CD38", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD38.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "PTPRC")
p2 <- FeaturePlot(seurat.object.subset, features="PTPRC", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_PTPRC.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD38_FPR2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "CD38")
p2 <- FeaturePlot(seurat.object.subset, features="CD38", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD38.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "FPR2")
p2 <- FeaturePlot(seurat.object.subset, features="FPR2", slot="scale.data")
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
Idents(seurat.object.subset) <- "Type"

# SARS.genes <- c("ORF1AB-1-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
p1 <- plot_density(seurat.object.subset, "CD4")
p2 <- FeaturePlot(seurat.object.subset, "CD8A")
p5 <- CombinePlots(plots = list(p1, p2))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("UMAP_KernelDensity_CD4-CD8A.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots
p3 <- plot_density(seurat.object.subset, c("CD4", "CD8A"))
p5 <- p3 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_CD4-CD8A_Joint.pdf", plot = p5, device = "pdf", width = 5, height = 7, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "CD8A"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_major.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("TRGV9", "TRDV2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_gd.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("TRAV10", "TRAJ18"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MAIT.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD3G", "FCGR3A","NCAM1", "NCR1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_NK.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")

	#	Error in linbin2D.ks(x, gpoints[[1]], gpoints[[2]], w = w) : 
	#	  NA/NaN/Inf in foreign function call (arg 10)
	#	In addition: Warning message:
	#	In ks.defaults(x = x, w = w, binned = binned, bgridsize = bgridsize,  :
	#	  Weights don't sum to sample size - they have been scaled accordingly


## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("IFNG", "TBX21", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH1.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 16, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("GATA3", "IL4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH2.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("RORC", "IL17A", "IL17F", "IL21"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH17.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")




## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(seurat.object.subset, features = lymphoid.lineage)
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
feature.plot <- DotPlot(seurat.object.subset, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(seurat.object.subset, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

macrophage.lineage <- c("CD14", "CD68", "CSF1R", "ADGRE1", "SIGLEC1", "CD38", "GPR18", "FPR2", "MYC", "EGR2")
feature.plot <- DotPlot(seurat.object.subset, features = macrophage.lineage)
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

feature.plot <- DotPlot(seurat.object.subset, features = macrophage.lineage.markers.2)
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

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=t.lineages)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, slot="counts", assay="Spatial")
ggsave("t.lineages_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, assay="Spatial")
ggsave("t.lineages_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=macrophage.lineage.markers.2)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("macrophage.lineage.markers.2_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, assay="Spatial")
ggsave("macrophage.lineage.markers.2_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=lymphoid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("lymphoid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, assay="Spatial")
ggsave("lymphoid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=myeloid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("myeloid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, assay="Spatial")
ggsave("myeloid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
# seurat.object.subset <- DietSeurat(seurat.object.subset, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

############## Activation gating strategy https://www.nature.com/articles/s41590-021-01049-2/figures/10
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD3D", "NCAM1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TorNK_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "KLRB1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRB1-CD161.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "CD38"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD38.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "CD69"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD69.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD4", "KLRK1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRK1-NKG2D.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD4", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("FOXP3", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-FOXP3-CD4.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("KIT", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("TNF", "TNFRSF8"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.TNFRSF8.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")


p4 <- plot_density(seurat.object.subset, c("CD14", "ITGAM"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("SIGLEC1", "FOLR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("SIGLEC1", "DDX3Y"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC-fetalsex_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
	## Placenta-associated maternal macrophages PAMMS (FOLR2 neg, HLA positive); HBCs FOLR2+HLA- ; Thomas et al. Naomi McGovern, JEM 2020
p4 <- plot_density(seurat.object.subset, c("FOLR2", "HLA-DRA"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMMs_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM1_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM2_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("ITGAM", "CD14"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M0.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD80", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD11B", "MYC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD163", "MRC1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2a.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("HLA-DRA", "IL6"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD14", "CD163"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2c.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD14", "CD1C"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD14", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p1 <- plot_density(seurat.object.subset, "THBD")
p2 <- FeaturePlot(seurat.object.subset, "THBD")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "ITGB3")
p2 <- FeaturePlot(seurat.object.subset, "ITGB3")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_Megakaryocyte.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "CD1C")
p2 <- FeaturePlot(seurat.object.subset, "CD1C")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "IL6")
p2 <- FeaturePlot(seurat.object.subset, "IL6")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "TNF")
p2 <- FeaturePlot(seurat.object.subset, "TNF")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "MYC")
p2 <- FeaturePlot(seurat.object.subset, "MYC")
p5 <- CombinePlots(plots = list(p1, p2))


ggsave("UMAP_M2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "FOLR2")
p2 <- FeaturePlot(seurat.object.subset, "FOLR2")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_HBCs-FOLR2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "HLA-DRA")
p2 <- FeaturePlot(seurat.object.subset, "HLA-DRA")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_PAMM-HLA-DRA.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "CD33")
p2 <- FeaturePlot(seurat.object.subset, "CD33")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeliodProgenitor-M0.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p4 <- plot_density(seurat.object.subset, c("KIT", "CD33"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_MyeliodProgenitor-M0-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

setwd("..")


############################################################# ############################################################# #############################################################
############################################################# Subset analysis of NegativeControl #############################################################
############################################################# ############################################################# #############################################################


dir.create("Subset-NegativeControl")
setwd("Subset-NegativeControl")
Idents(seurat.object)
Idents(seurat.object) <- "AnalysisCohort"

seurat.object.subset <- subset(x = seurat.object, idents = c("NegativeControl"))


p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "Subtype", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Subtype')
ggsave("UMAP-Subtype.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p4 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "SampleCode", cols = mypal3, raster=FALSE) + labs(title = NULL, color='SampleCode')
ggsave("UMAP-SampleCode.pdf".subset, plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")
p1 <- DimPlot(seurat.object.subset, reduction = "umap", group.by = "AnalysisCohort")
ggsave("UMAP-AnalysisCohort.pdf", plot = p1, device = "pdf", width = 7, height = 5, units = "in")

Idents(seurat.object.subset) <- "Type"
de_markers <- FindAllMarkers(seurat.object.subset, features = all.genes, assay = "Spatial", only.pos = FALSE, min.pct = 0.01, logfc.threshold =  2)
write.table(de_markers, "DEGs_byType_2log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, raster = FALSE)
ggsave("top5_markers_logfc_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")
top5.heatmap <- DoHeatmap(seurat.object.subset, features = top5$gene, slot="counts", raster = FALSE)
ggsave("top5_markers_counts_heatmap.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 9, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top10.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top15.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object.subset, features = unique.top2)
png(file=paste0("top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Type")
dev.off()
Idents(seurat.object.subset) <- "Type"
## use RNA for all DE analysis and plots
# seurat.object.subset <- ScaleData(seurat.object.subset, features = all.genes, assay = "SCT")
# seurat.object.subset <- ScaleData(seurat.object.subset, features = all.genes, assay = "Spatial")
DefaultAssay(seurat.object.subset) <- "Spatial"
Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("ITGAM", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_ITGAM_THBD.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "ITGAM")
p2 <- FeaturePlot(seurat.object.subset, features="ITGAM", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "THBD")
p2 <- FeaturePlot(seurat.object.subset, features="THBD", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_THBD.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "IFNG"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD4_IFNG.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "CD4")
p2 <- FeaturePlot(seurat.object.subset, features="CD4", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD4.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "IFNG")
p2 <- FeaturePlot(seurat.object.subset, features="IFNG", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_IFNG.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD38", "PTPRC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD38_PTPRC.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "CD38")
p2 <- FeaturePlot(seurat.object.subset, features="CD38", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD38.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "PTPRC")
p2 <- FeaturePlot(seurat.object.subset, features="PTPRC", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_PTPRC.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

Idents(seurat.object.subset) <- "Type"
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("all_CD38_FPR2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p1 <- plot_density(seurat.object.subset, "CD38")
p2 <- FeaturePlot(seurat.object.subset, features="CD38", slot="scale.data")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("all_CD38.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

p1 <- plot_density(seurat.object.subset, "FPR2")
p2 <- FeaturePlot(seurat.object.subset, features="FPR2", slot="scale.data")
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
Idents(seurat.object.subset) <- "Type"

# SARS.genes <- c("ORF1AB-1-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
p1 <- plot_density(seurat.object.subset, "CD4")
p2 <- FeaturePlot(seurat.object.subset, "CD8A")
p5 <- CombinePlots(plots = list(p1, p2))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("UMAP_KernelDensity_CD4-CD8A.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots
p3 <- plot_density(seurat.object.subset, c("CD4", "CD8A"))
p5 <- p3 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_CD4-CD8A_Joint.pdf", plot = p5, device = "pdf", width = 5, height = 7, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "CD8A"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_major.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("TRGV9", "TRDV2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_gd.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("TRAV10", "TRAJ18"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MAIT.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD3G", "FCGR3A","NCAM1", "NCR1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_NK.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")

	#	Error in linbin2D.ks(x, gpoints[[1]], gpoints[[2]], w = w) : 
	#	  NA/NaN/Inf in foreign function call (arg 10)
	#	In addition: Warning message:
	#	In ks.defaults(x = x, w = w, binned = binned, bgridsize = bgridsize,  :
	#	  Weights don't sum to sample size - they have been scaled accordingly


## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("IFNG", "TBX21", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH1.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 16, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("GATA3", "IL4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH2.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("RORC", "IL17A", "IL17F", "IL21"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH17.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")




## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(seurat.object.subset, features = lymphoid.lineage)
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
feature.plot <- DotPlot(seurat.object.subset, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(seurat.object.subset, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

macrophage.lineage <- c("CD14", "CD68", "CSF1R", "ADGRE1", "SIGLEC1", "CD38", "GPR18", "FPR2", "MYC", "EGR2")
feature.plot <- DotPlot(seurat.object.subset, features = macrophage.lineage)
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

feature.plot <- DotPlot(seurat.object.subset, features = macrophage.lineage.markers.2)
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

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=t.lineages)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, slot="counts", assay="Spatial")
ggsave("t.lineages_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, assay="Spatial")
ggsave("t.lineages_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=macrophage.lineage.markers.2)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, slot="counts", assay="Spatial")
ggsave("macrophage.lineage.markers.2_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, assay="Spatial")
ggsave("macrophage.lineage.markers.2_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=lymphoid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("lymphoid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, assay="Spatial")
ggsave("lymphoid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object.subset, assays="Spatial", slot="counts", return.seurat = TRUE, features=myeloid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, slot="counts", assay="Spatial")
ggsave("myeloid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, assay="Spatial")
ggsave("myeloid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
# seurat.object.subset <- DietSeurat(seurat.object.subset, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

############## Activation gating strategy https://www.nature.com/articles/s41590-021-01049-2/figures/10
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD3D", "NCAM1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TorNK_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "KLRB1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRB1-CD161.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "CD38"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD38.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(seurat.object.subset, c("CD4", "CD69"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD69.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD4", "KLRK1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRK1-NKG2D.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD4", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("FOXP3", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-FOXP3-CD4.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("KIT", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("TNF", "TNFRSF8"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.TNFRSF8.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")


p4 <- plot_density(seurat.object.subset, c("CD14", "ITGAM"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("SIGLEC1", "FOLR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("SIGLEC1", "DDX3Y"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC-fetalsex_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
	## Placenta-associated maternal macrophages PAMMS (FOLR2 neg, HLA positive); HBCs FOLR2+HLA- ; Thomas et al. Naomi McGovern, JEM 2020
p4 <- plot_density(seurat.object.subset, c("FOLR2", "HLA-DRA"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMMs_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM1_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM2_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("ITGAM", "CD14"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M0.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD80", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD11B", "MYC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD163", "MRC1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2a.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("HLA-DRA", "IL6"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD14", "CD163"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2c.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD14", "CD1C"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(seurat.object.subset, c("CD14", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p1 <- plot_density(seurat.object.subset, "THBD")
p2 <- FeaturePlot(seurat.object.subset, "THBD")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "ITGB3")
p2 <- FeaturePlot(seurat.object.subset, "ITGB3")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_Megakaryocyte.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "CD1C")
p2 <- FeaturePlot(seurat.object.subset, "CD1C")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "IL6")
p2 <- FeaturePlot(seurat.object.subset, "IL6")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "TNF")
p2 <- FeaturePlot(seurat.object.subset, "TNF")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "MYC")
p2 <- FeaturePlot(seurat.object.subset, "MYC")
p5 <- CombinePlots(plots = list(p1, p2))


ggsave("UMAP_M2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "FOLR2")
p2 <- FeaturePlot(seurat.object.subset, "FOLR2")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_HBCs-FOLR2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "HLA-DRA")
p2 <- FeaturePlot(seurat.object.subset, "HLA-DRA")
p5 <- CombinePlots(plots = list(p1, p2))

ggsave("UMAP_PAMM-HLA-DRA.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(seurat.object.subset, "CD33")
p2 <- FeaturePlot(seurat.object.subset, "CD33")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeliodProgenitor-M0.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p4 <- plot_density(seurat.object.subset, c("KIT", "CD33"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_MyeliodProgenitor-M0-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

setwd("..")


############################################################# ############################################################# #############################################################
############################################################# Make SpatialPlots from Top Markers identified above #############################################################
############################################################# ############################################################# #############################################################

dir.create("TopFeature-SpatialPlots")
setwd("TopFeature-SpatialPlots")
Idents(seurat.object)
Idents(seurat.object) <- "AnalysisCohort"

feature.list <- c("PSG7", "CSH2", "GDF15", "KISS1")

image.list <- c("S01", "S03.1", "S04.2", 
    "S24.13", "S25.14", "S26.15",
 "S15.3", "S20.8", "S21.9", 
 "S16.4", "S17.5", "S18.6", 
 "S19.7", "S22.10", "S23a.11", "S23b.12")

#p5 <- SpatialFeaturePlot(seurat.object, features = feature.list, images=image.list, ncol=4)
#ggsave("all_SpatialPlots-cropped-NDSP-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-NDSP-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S15.3", "S20.8", "S21.9")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NotDetected_SpatialPlots-cropped-NDSP-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S16.4", "S17.5", "S18.6", "S19.7", "S22.10")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("SparsePositive_SpatialPlots-cropped-NDSP-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S23a.11", "S23b.12")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=2)
ggsave("HighPositive_SpatialPlots-cropped-NDSP-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

feature.list <- c("IFITM3", "IL32", "IL1RN", "CXCL8", "CCL20", "IL24", "ISG15", "IFIT3", "IL1RL1", "CXCL10", "MT2A", "TIMP1", "GBP1", "S100A9")
#p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=4)
#ggsave("all_SpatialPlots-cropped-HP-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-HP-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S15.3", "S20.8", "S21.9")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NotDetected_SpatialPlots-cropped-HP-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S16.4", "S17.5", "S18.6", "S19.7", "S22.10")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("SparsePositive_SpatialPlots-cropped-HP-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("S23a.11", "S23b.12")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=2)
ggsave("HighPositive_SpatialPlots-cropped-HP-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")



feature.list <- c("TIMP1")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-TIMP1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("S15.3", "S20.8", "S21.9")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NotDetected_SpatialPlots-cropped-TIMP1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("S16.4", "S17.5", "S18.6", "S19.7", "S22.10")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("SparsePositive_SpatialPlots-cropped-TIMP1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("S23a.11", "S23b.12")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=2)
ggsave("HighPositive_SpatialPlots-cropped-TIMP1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

Idents(seurat.object) <- "AnalysisCohort"
library.averages.heatmap <- VlnPlot(seurat.object, features = "TIMP1", assay="Spatial", group.by="AnalysisCohort")
library.averages.heatmap
ggsave("vlnplot_TIMP1_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")


Idents(seurat.object) <- "AnalysisCohort"
library.averages.heatmap <- VlnPlot(seurat.object, features = "KISS1", assay="Spatial", group.by="AnalysisCohort")
library.averages.heatmap
ggsave("vlnplot_KISS1_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

feature.list <- c("KISS1")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-KISS1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("S15.3", "S20.8", "S21.9")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NotDetected_SpatialPlots-cropped-KISS1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("S16.4", "S17.5", "S18.6", "S19.7", "S22.10")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("SparsePositive_SpatialPlots-cropped-KISS1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("S23a.11", "S23b.12")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=2)
ggsave("HighPositive_SpatialPlots-cropped-KISS1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


feature.list <- c("PSG7", "CSH2", "GDF15", "KISS1")
library.averages.heatmap <- VlnPlot(seurat.object, features = feature.list, assay="Spatial", group.by="AnalysisCohort")
library.averages.heatmap
ggsave("vlnplot_NotDetected_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")


feature.list <- c("IFITM3", "IL32", "IL1RN", "CXCL8", "CCL20", "IL24", "ISG15", "IFIT3", "IL1RL1", "CXCL10", "MT2A", "TIMP1", "GBP1", "S100A9")
library.averages.heatmap <- VlnPlot(seurat.object, features = feature.list, assay="Spatial", group.by="AnalysisCohort")
library.averages.heatmap
ggsave("vlnplot_Cytokines_orig.ident.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")


seurat.object.subset <- subset(x = seurat.object.subset, idents = c("NegativeControl"))
Idents(seurat.object.subset) <- "SampleCodes"
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object.subset, features = "feature.list", images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-topmarkers_NegativeControl_subset.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

seurat.object.subset <- subset(x = seurat.object.subset, idents = c("NotDetected"))
Idents(seurat.object.subset) <- "SampleCodes"
image.list <- c("S15.3", "S20.8", "S21.9")
p5 <- SpatialPlot(seurat.object.subset, features = "feature.list", images=image.list, ncol=3)
ggsave("NotDetected_SpatialPlots-cropped-topmarkers_NotDetected_subset.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

seurat.object.subset <- subset(x = seurat.object.subset, idents = c("SparsePositive"))
Idents(seurat.object.subset) <- "SampleCodes"
image.list <- c("S16.4", "S17.5", "S18.6", "S19.7", "S22.10")
p5 <- SpatialPlot(seurat.object.subset, features = "feature.list", images=image.list, ncol=3)
ggsave("SparsePositive_SpatialPlots-cropped-topmarkers_SparsePositive_subset.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

seurat.object.subset <- subset(x = seurat.object.subset, idents = c("HighPositive"))
Idents(seurat.object.subset) <- "SampleCodes"
image.list <- c("S01", "S03.1", "S04.2", 
    "S24.13", "S25.14", "S26.15",
 "S15.3", "S20.8", "S21.9", 
 "S16.4", "S17.5", "S18.6", 
 "S19.7", "S22.10", "S23a.11", "S23b.12")
image.list <- c("S23a.11", "S23b.12")
p5 <- SpatialPlot(seurat.object.subset, features = "feature.list", images=image.list, ncol=2)
ggsave("HighPositive_SpatialPlots-cropped-topmarkers_HighPositive_subset.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

############################################################# ############################################################# #############################################################
############################################################# Make SpatialPlots of macrophages #############################################################
############################################################# ############################################################# #############################################################

dir.create("Subset-macrophage")
setwd("Subset-macrophage")
Idents(seurat.object)
Idents(seurat.object) <- "Type"
levels(seurat.object)

seurat.object.subset <- subset(x = seurat.object, idents = c("Macrophage"))


