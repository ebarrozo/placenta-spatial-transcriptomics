## human_seurat_interated_v3.R

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

## integrating our visium spatial transcriptomics data with Tsang et al (2017) human placenta data (term placenta including samples with pre-eclampsia)
	## v1 was first 4 visium samples
	## v2 has all visium samples + tsang and saha
	## v3 is also adds Wayne-State dbGaP placental scRNA-seq data; removes saha data and keeps Preeclampsia tsang data; adds Yang-2021 dataset

## Analyze data on AagaardLab3
# 

set.seed(seed=1)
setwd("/home/ebarrozo/visium/results")
library(dplyr)
library(Seurat)	 # Seurat_4.0.3
library(sctransform)
library(ggplot2)
library(cowplot)
library(SeuratDisk)
library(dplyr)
library(patchwork)
# BiocManager::install("glmGamPoi")
library(glmGamPoi)
library(ggpubr)


###### load yang et al 2021 dataset
setwd("/home/ebarrozo/placenta-scrnaseq/yang/results/seurat-v1")
load("Yang_data-integrated_umap_annotated-v1.RData")
seurat.object.yang <- seurat.object
rm(seurat.object)

setwd("/home/ebarrozo/visium/results/seurat_human_v2")
load("/home/ebarrozo/visium/results/seurat_human_v2/annotated/human_spatial_data-integrated-annotated_v1.RData")
new.metadata <- c("PV", "BP", "CAM", "PVBP", 
	"PVBP", "PVBP", "PVBP", "PVBP", 
	"PVBP", "PVBP", "PVBP", "PVBP", 
	"PVBP", "PVBP", "PVBP", "PVBP")
Idents(seurat.object) <- "orig.ident"
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
## add cluster Subtypes as a metadata column
seurat.object$Site <- Idents(seurat.object)

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("Control", "Control", "Control", "SARS-CoV-2.ax", "SARS-CoV-2.sx", "SARS-CoV-2.ax", "SARS-CoV-2.ax", "SARS-CoV-2.ax", "SARS-CoV-2.sx", "SARS-CoV-2.sx", "SARS-CoV-2.sx", "SARS-CoV-2.sx", "SARS-CoV-2.sx", "Control","Control","Control")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Condition <- Idents(seurat.object)
Idents(seurat.object) <- "Condition"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("Barrozo-2022", "Barrozo-2022", "Barrozo-2022", "Barrozo-2022", "Barrozo-2022", "Barrozo-2022", "Barrozo-2022", "Barrozo-2022", "Barrozo-2022", "Barrozo-2022", "Barrozo-2022", "Barrozo-2022", "Barrozo-2022", "Barrozo-2022","Barrozo-2022","Barrozo-2022")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$etal <- Idents(seurat.object)
Idents(seurat.object) <- "etal"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("Spatial", "Spatial", "Spatial", "Spatial", "Spatial", "Spatial", "Spatial", "Spatial", "Spatial", "Spatial", "Spatial", "Spatial", "Spatial", "Spatial","Spatial","Spatial")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Platform <- Idents(seurat.object)
Idents(seurat.object) <- "Platform"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("F", "F", "F", "F", 
	"M", "F", "M", "M", 
	"F", "F", "M", "F", 
	"F", "M","M","M")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$FetalSex <- Idents(seurat.object)
Idents(seurat.object) <- "Platform"



	## see ## visium_human_seurat_v1.R
seurat.object.spatial <- seurat.object
rm(seurat.object)

setwd("/home/ebarrozo/placenta-scrnaseq/tsang-saha_merged")
load("/home/ebarrozo/placenta-scrnaseq/tsang-saha_merged/tsang-saha_merged_data-annotated.RData")
	## see  Tsang2017-EGAS00001002449-seurat_v3.R for QC pre-processing cluster annotations

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
[1] "PE1"         "PE2"         "PE3"         "PE4"         "PN1"         "PN2"        
 [7] "PN3C"        "PN3P"        "PN4C"        "PN4P"        "six.weeks"   "seven.weeks"
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE)

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
cluster.Type <- c("Preeclampsia", "Preeclampsia", "Preeclampsia", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control")
names(cluster.Type) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Type)
seurat.object$Condition <- Idents(seurat.object)
Idents(seurat.object) <- "Condition"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
cluster.Type <- c("PVBP", "PVBP", "PVBP", "PVBP", "PVBP", "PVBP", "PVBP", "PVBP", "PVBP", "PVBP", "PVBP", "PVBP")
names(cluster.Type) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Type)
seurat.object$Site <- Idents(seurat.object)
Idents(seurat.object) <- "Site"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
cluster.Type <- c("scRNA-seq", "scRNA-seq", "scRNA-seq", "scRNA-seq", "scRNA-seq", "scRNA-seq", "scRNA-seq", "scRNA-seq", "scRNA-seq", "scRNA-seq", "scRNA-seq", "scRNA-seq")
names(cluster.Type) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Type)
seurat.object$Platform <- Idents(seurat.object)
Idents(seurat.object) <- "Platform"

ncol(seurat.object)
	# 41505
DefaultAssay(seurat.object) <- "SCT"
Idents(seurat.object) <- "source"
levels(seurat.object)

seurat.object <- subset(x = seurat.object, idents = c("tsang"))
ncol(seurat.object)
	# 32248
Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
cluster.Type <- c("Tsang-2017", "Tsang-2017", "Tsang-2017","Tsang-2017","Tsang-2017","Tsang-2017", "Tsang-2017", "Tsang-2017","Tsang-2017","Tsang-2017")
names(cluster.Type) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.Type)
seurat.object$etal <- Idents(seurat.object)
Idents(seurat.object) <- "etal"

seurat.object.tsang <- seurat.object
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

setwd("/home/ebarrozo/visium/results/seurat_human_v2")
dir.create("integration-v3")
setwd("integration-v3")
setwd("/home/ebarrozo/visium/results/seurat_human_v2/integration-v3")

############################################################# ############################################################# #############################################################
############################################################# Integration: merge all data using reciprical PCA anchors #############################################################
############################################################# ############################################################# #############################################################
## Merge objects
all_merged <- merge(x = seurat.object, y = c(seurat.object.spatial, seurat.object.tsang, seurat.object.yang), merge.data = TRUE, project = "integrated")
rm(seurat.object.spatial)
rm(seurat.object)
rm(seurat.object.tsang)
rm(seurat.object.yang)

ncol(all_merged)
	# 291,871 transcriptomes # 49 GB
		## ncol minus 17927 spatial transcriptomes = 273,944
Idents(all_merged) <- "orig.ident"
t1<-table(Idents(all_merged))
write.table(t1, "orig.ident.counts.txt", sep="\t")
Idents(all_merged) <- "Type"
t1<-table(Idents(all_merged))
write.table(t1, "Type.counts.txt", sep="\t")		
Idents(all_merged) <- "Platform"
t1<-table(Idents(all_merged))
write.table(t1, "Platform.counts.txt", sep="\t")
Idents(all_merged) <- "Site"
t1<-table(Idents(all_merged))
write.table(t1, "Site.counts.txt", sep="\t")
Idents(all_merged) <- "Condition"
t1<-table(Idents(all_merged))
write.table(t1, "Condition.counts.txt", sep="\t")
Idents(all_merged) <- "etal"
t1<-table(Idents(all_merged))
write.table(t1, "etal.counts.txt", sep="\t")				
# Pique-Regi-2019    Pique-Regi-2020    Garcia-Flores-2022       Barrozo-2022         Tsang-2017 	 Yang-2021 
#        50,710               5,232             147,906              17,927              32,248 		37,848

all.genes <- rownames(all_merged@assays[["RNA"]])
	## 53k genes? includes ncRNAs like C19ORF70, "AP001610.1"      "H2BFS.1"         "AP001065.1"  etc
# all.genes.spatial <- rownames(all_merged@assays[["Spatial"]])
	## 36k genes
			## use .spatial as all.genes
#all.genes <- union(all.genes, all.genes.spatial)
all.genes <- rownames(all_merged@assays[["Spatial"]])
options(future.globals.maxSize = 250000000000)

rm(object.11)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(umap.combined)
rm(all.genes.spatial)
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


setwd("/home/ebarrozo/visium/results/seurat_human_v2/integration-v3")

# 	save.image("human-merged_v3.RData")
load("human-merged_v3.RData")



all_merged <- PercentageFeatureSet(all_merged, features = Spike.genes, col.name = "percent.Spike", assay="SCT")
all_merged <- PercentageFeatureSet(all_merged, features = HB.genes, col.name = "percent.Hemo", assay="SCT")
all_merged <- PercentageFeatureSet(all_merged, features = "DDX3Y", col.name = "percent.DDX3Y", assay="SCT")

SARS.receptors <- c("ACE2", "TMPRSS2", "CTSB", "CTSL", "HBA1", "HBB", "HBD", "HBZ", "HBQ1", "ICAM4", "DDX3Y", "SPP1")
SARS.genes.receptors.hb.fetalsex.mac <- union(SARS.genes, SARS.receptors)


Idents(all_merged) <- "orig.ident"
new.metadata <- c("M", "M", "F", "M", "F", "F", "F", "F", "F", "M", "M", "M", "M", "M", "F", "F", "M", "F", "F", "F", "M", "M", "M", "M", "M", "M", "M", "M", "M", "F", "F", "F", "M", "F", "M", "F", "F", "M", "M", "F", "M", "M", "M", "M", "M", "F", "F", "F", "F", "M", "F", "M", "M", "F", "F", "M", "F", "M", "M", "M", "M", "M", "M", "F", "F", "M", "M", "F", "F", "M", "F", "F", "F", "F", "F", "F", "F", "F", "F")
names(new.metadata) <- levels(all_merged)
all_merged <- RenameIdents(all_merged, new.metadata)
all_merged$FetalSex <- Idents(all_merged)







all_merged <- SCTransform(all_merged, method = "glmGamPoi", assay = "SCT", return.only.var.genes=TRUE, conserve.memory = TRUE, verbose = TRUE)
# all_merged <- SCTransform(all_merged, method = "glmGamPoi", assay = "SCT", return.only.var.genes=FALSE, verbose = T)

	## Keeps crashing here, so using new method (glmGamPoi) and conserve.memory option.
			## 49gb object; 103 gb rstudio session

all_merged <- PercentageFeatureSet(all_merged, features = SARS.genes, col.name = "percent.viral", assay="SCT")
	# error
all_merged <- PercentageFeatureSet(all_merged, features = Spike.genes, col.name = "percent.Spike", assay="SCT")
all_merged <- PercentageFeatureSet(all_merged, features = HB.genes, col.name = "percent.Hemo", assay="SCT")
all_merged <- PercentageFeatureSet(all_merged, features = "DDX3Y", col.name = "percent.DDX3Y", assay="SCT")

SARS.receptors <- c("ACE2", "TMPRSS2", "CTSB", "CTSL", "HBA1", "HBB", "HBD", "HBZ", "HBQ1", "ICAM4", "DDX3Y", "SPP1")
SARS.genes.receptors.hb.fetalsex.mac <- union(SARS.genes, SARS.receptors)

Idents(all_merged) <- "orig.ident"
levels(all_merged)
 [1] "SRR10166475" "SRR10166477" "SRR10166480" "SRR10166483" "SRR10166486" "SRR10166489"
 [7] "SRR10166492" "SRR10166495" "SRR10166498" "SRR12192596" "SRR12192597" "SRR12192598"
[13] "SRR12192599" "SRR10166484" "SRR10166481" "SRR10166496" "SRR10166478" "SRR10166487"
[19] "SRR10166490" "SRR10166493" "SRR17065313" "SRR17065314" "SRR17065315" "SRR17065316"
[25] "SRR17065317" "SRR17065318" "SRR17065319" "SRR17065320" "SRR10166479" "SRR10166488"
[31] "SRR10166494" "SRR10166499" "SRR10166476" "SRR10166482" "SRR10166485" "SRR10166491"
[37] "SRR10166497" "SRR17065305" "SRR17065306" "SRR17065307" "SRR17065308" "SRR17065309"
[43] "SRR17065310" "SRR17065311" "SRR17065312" "S01"         "S03"         "S04"        
[49] "S15"         "S16"         "S17"         "S18"         "S19"         "S20"        
[55] "S21"         "S22"         "S23a"        "S23b"        "S24"         "S25"        
[61] "S26"         "PE1"         "PE2"         "PE3"         "PE4"         "PN1"        
[67] "PN2"         "PN3C"        "PN3P"        "PN4C"        "PN4P"        "C1"         
[73] "C2"          "C3"          "C4"          "P1"          "P2"          "G1"         
[79] "G2"  




Idents(all_merged) <- "etal"
cluster.averages <- AverageExpression(all_merged, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes.receptors.hb.fetalsex.mac)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, slot="counts", assay="SCT", size=4, angle = 45) + NoLegend()
	#  The following features were omitted as they were not found in the scale.data slot for the SCT assay: ORF7B, ORF6, E
ggsave("SARS.genes.receptors.hb.fetalsex.mac_mc_heatmap_etal.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, assay="SCT", size=4, angle = 45) + NoLegend()
ggsave("SARS.genes.receptors.hb.fetalsex.mac_FC_heatmap_etal.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
Idents(all_merged) <- "orig.ident"
cluster.averages <- AverageExpression(all_merged, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes.receptors.hb.fetalsex.mac)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, slot="counts", assay="SCT", size=2, angle = 90) + NoLegend()
ggsave("SARS.genes.receptors.hb.fetalsex.mac_mc_heatmap_orig.ident.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, assay="SCT", size=2, angle = 90) + NoLegend()
ggsave("SARS.genes.receptors.hb.fetalsex.mac_FC_heatmap_orig.ident.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
Idents(all_merged) <- "orig.ident"
sex.genes <- c("SRY", "DDX3Y")
cluster.averages <- AverageExpression(all_merged, assays="SCT", slot="counts", return.seurat = TRUE, features=sex.genes)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = sex.genes, raster = FALSE, slot="counts", assay="SCT", size=3, angle = 90) + NoLegend()
ggsave("fetalsex_mc_heatmap_orig.ident.png", plot = library.averages.heatmap, device = "png", width = 11, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = sex.genes, raster = FALSE, assay="SCT", size=3, angle = 90) + NoLegend()
ggsave("fetalsex_FC_heatmap_orig.ident.png", plot = library.averages.heatmap, device = "png", width = 11, height = 5, units = "in")

write.csv(cluster.averages@assays[["SCT"]]@counts, file = "orig.ident_DDX3Y_meancounts.csv")

Idents(all_merged) <- "orig.ident"
new.metadata <- c("M", "M", "F", "M", "F", "F", "F", "F", "F", "M", "M", "M", "M", "M", "F", "F", "M", "F", "F", "F", "M", "M", "M", "M", "M", "M", "M", "M", "M", "F", "F", "F", "M", "F", "M", "F", "F", "M", "M", "F", "M", "M", "M", "M", "M", "F", "F", "F", "F", "M", "F", "M", "M", "F", "F", "M", "F", "M", "M", "M", "M", "M", "M", "F", "F", "M", "M", "F", "F", "M", "F", "F", "F", "F", "F", "F", "F", "F", "F")
names(new.metadata) <- levels(all_merged)
all_merged <- RenameIdents(all_merged, new.metadata)
all_merged$FetalSex <- Idents(all_merged)

Idents(all_merged) <- "Condition"
cluster.averages <- AverageExpression(all_merged, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes.receptors.hb.fetalsex.mac)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, slot="counts", assay="SCT", size=4, angle = 45) + NoLegend()
ggsave("SARS.genes.receptors.hb.fetalsex.mac_mc_heatmap_condition.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, assay="SCT", size=4, angle = 45) + NoLegend()
ggsave("SARS.genes.receptors.hb.fetalsex.mac_FC_heatmap_condition.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
Idents(all_merged) <- "Site"
cluster.averages <- AverageExpression(all_merged, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes.receptors.hb.fetalsex.mac)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, slot="counts", assay="SCT", size=4, angle = 45) + NoLegend()
ggsave("SARS.genes.receptors.hb.fetalsex.mac_mc_heatmap_site.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, assay="SCT", size=4, angle = 45) + NoLegend()
ggsave("SARS.genes.receptors.hb.fetalsex.mac_FC_heatmap_site.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
Idents(all_merged) <- "Platform"
cluster.averages <- AverageExpression(all_merged, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes.receptors.hb.fetalsex.mac)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, slot="counts", assay="SCT", size=4, angle = 45) + NoLegend()
ggsave("SARS.genes.receptors.hb.fetalsex.mac_mc_heatmap_Platform.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, assay="SCT", size=4, angle = 45) + NoLegend()
ggsave("SARS.genes.receptors.hb.fetalsex.mac_FC_heatmap_Platform.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
Idents(all_merged) <- "Type"
cluster.averages <- AverageExpression(all_merged, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes.receptors.hb.fetalsex.mac)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, slot="counts", assay="SCT", size=4, angle = 45) + NoLegend()
ggsave("SARS.genes.receptors.hb.fetalsex.mac_mc_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 14, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, assay="SCT", size=4, angle = 45) + NoLegend()
ggsave("SARS.genes.receptors.hb.fetalsex.mac_FC_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 14, height = 8, units = "in")
rm(library.averages.heatmap)
rm(cluster.averages)


# all_merged <- DietSeurat(all_merged, counts=TRUE, data=TRUE, scale.data = FALSE)
		## NEEEEED scale data for PCA below
## from 48gb to 41

Idents(all_merged) <- "etal"
levels(all_merged)
# save.image("human-var.genes-scaled_v3.RData")



set.seed(seed=1)
setwd("/home/ebarrozo/visium/results/seurat_human_v2/integration-v3")
library(dplyr)
library(Seurat)	 # Seurat_4.0.3
library(sctransform)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(glmGamPoi)
library(ggpubr)

# load("human-var.genes-scaled_v3.RData")


setwd("/home/ebarrozo/visium/results/seurat_human_v2/integration-v3")

# 	save.image("human-merged_v3.RData")
# # 	load("human-merged_v3.RData")



all_merged <- PercentageFeatureSet(all_merged, features = Spike.genes, col.name = "percent.Spike", assay="SCT")
all_merged <- PercentageFeatureSet(all_merged, features = HB.genes, col.name = "percent.Hemo", assay="SCT")
all_merged <- PercentageFeatureSet(all_merged, features = "DDX3Y", col.name = "percent.DDX3Y", assay="SCT")

SARS.receptors <- c("ACE2", "TMPRSS2", "CTSB", "CTSL", "HBA1", "HBB", "HBD", "HBZ", "HBQ1", "ICAM4", "DDX3Y", "SPP1")
SARS.genes.receptors.hb.fetalsex.mac <- union(SARS.genes, SARS.receptors)


Idents(all_merged) <- "orig.ident"
new.metadata <- c("M", "M", "F", "M", "F", "F", "F", "F", "F", "M", "M", "M", "M", "M", "F", "F", "M", "F", "F", "F", "M", "M", "M", "M", "M", "M", "M", "M", "M", "F", "F", "F", "M", "F", "M", "F", "F", "M", "M", "F", "M", "M", "M", "M", "M", "F", "F", "F", "F", "M", "F", "M", "M", "F", "F", "M", "F", "M", "M", "M", "M", "M", "M", "F", "F", "M", "M", "F", "F", "M", "F", "F", "F", "F", "F", "F", "F", "F", "F")
names(new.metadata) <- levels(all_merged)
all_merged <- RenameIdents(all_merged, new.metadata)
all_merged$FetalSex <- Idents(all_merged)








# [1] "Pique-Regi-2019"    "Pique-Regi-2020"    "Garcia-Flores-2022" "Barrozo-2022"      
# [5] "Tsang-2017"         "Yang-2021" 

# all_merged.list <- SplitObject(all_merged, split.by = "etal")

# Make a reference list
# reference.list <- all_merged.list[c("Tsang-2017", "Pique-Regi-2019", "Pique-Regi-2020", "Yang-2021", "Garcia-Flores-2022", "Barrozo-2022")]



# Select transcripts used for data integration
# reference.features <- SelectIntegrationFeatures(object.list = reference.list, nfeatures = 3000)
# write.table(reference.features, "all_merged.reference.features.txt", sep="\t")

# reference.list <- PrepSCTIntegration(object.list = reference.list, anchor.features = reference.features, verbose = T)

# Find integration anchors, this takes about 15-20 mins
# all_merged.anchors <- FindIntegrationAnchors(object.list = reference.list, normalization.method = "SCT", anchor.features = reference.features, verbose = T)


# seurat.object <- IntegrateData(anchorset = all_merged.anchors, normalization.method = "SCT")
# rm("all_merged.anchors")


## If CCA proves prohibitevly computationally expensive, utilize reciprical PCA
# all_merged<- ScaleData(all_merged, features = reference.features, assay = "SCT", verbose = TRUE)

dir.create("wo.integration.var.reg")
setwd("wo.integration.var.reg")
# all_merged<- ScaleData(all_merged, features = var.genes, assay = "SCT", verbose = TRUE)
all_merged <- SCTransform(all_merged, method = "glmGamPoi", assay = "SCT", return.only.var.genes=TRUE, do.scale=FALSE, conserve.memory = TRUE, vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"), verbose = TRUE)
# all_merged <- SCTransform(all_merged, method = "glmGamPoi", assay = "SCT", return.only.var.genes=TRUE, do.scale=TRUE, conserve.memory = TRUE, vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"), verbose = TRUE)

var.genes <- all_merged@assays[["SCT"]]@var.features
all_merged<- ScaleData(all_merged, features = var.genes, assay = "SCT", vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"), verbose = TRUE)

all_merged<- (all_merged, features = var.genes, assay = "SCT", slot = "scale.data", verbose = TRUE)

all_merged <- FindNeighbors(all_merged, features = "var.genes", dims = 1:30)
all_merged <- FindClusters(all_merged, resolution = 0.6)
all_merged <- RunUMAP(all_merged, dims=1:30)
	## 41 clusters with res=0.6; cell cycle effects; try again w/ regressing out percent.mt, percent.ribo, cc

# all_merged <- DietSeurat(all_merged, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

Idents(all_merged) <- "seurat_clusters"
t1<-table(Idents(all_merged))
write.table(t1, "seurat_clusters.counts.txt", sep="\t")
rm(t1)
# save.image("human-notintegrated-umap-scaled_v3.RData")

library(ggpubr)
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

p1 <- DimPlot(all_merged, reduction = "umap", group.by = "etal", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='et al.')
p1
ggsave("UMAP_XL-etal.pdf", plot = p1, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(all_merged, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(all_merged, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Profile')
p2
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 12, height = 6, units = "in")
p3 <- DimPlot(all_merged, reduction = "umap", group.by = "Site", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='Site')
ggsave("UMAP-Site.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p3 <- DimPlot(all_merged, reduction = "umap", group.by = "Platform", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='Platform')
ggsave("UMAP-Platform.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p4 <- DimPlot(all_merged, reduction = "umap", group.by = "Condition", cols = mypal3, raster=FALSE) + labs(title = NULL, color='Condition')
ggsave("UMAP-Condition.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")

umap.combined <- p1 + p2 + p3 + p4
ggsave("integrated_UMAP-FINAL.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("UMAP_XL.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")

p4 <- DimPlot(all_merged, reduction = "umap", group.by = "FetalSex", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='FetalSex')
ggsave("UMAP-FetalSex.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")
p4 <- FeaturePlot(all_merged, features="DDX3Y", reduction = "umap", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='DDX3Y')
ggsave("UMAP-DDX3Y.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")

#examine UMAPs with qc metrics 
p5 <- FeaturePlot(all_merged, features = 'nFeature_RNA')
p6 <- FeaturePlot(all_merged, features = 'nCount_RNA')
p7 <- FeaturePlot(all_merged, features = 'percent.mt')
p8 <- FeaturePlot(all_merged, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots-RNA.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")
#examine UMAPs with qc metrics 
p5 <- FeaturePlot(all_merged, features = 'nFeature_Spatial'
p6 <- FeaturePlot(all_merged, features = 'nCount_Spatial')
p7 <- FeaturePlot(all_merged, features = 'percent.mt')
p8 <- FeaturePlot(all_merged, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots-Spatial.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")
#examine UMAPs with qc metrics 
p5 <- DimPlot(all_merged, reduction = "umap", group.by = "Phase", cols = mypal3, raster=FALSE) + labs(title = NULL, color='Cell Cycle Phase')
p6 <- FeaturePlot(all_merged, features = 'percent.Hemo')
p7 <- FeaturePlot(all_merged, features = 'percent.Spike')
p8 <- FeaturePlot(all_merged, features = 'percent.viral')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots2.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")

Idents(all_merged) <- "etal"
# all_merged <- PercentageFeatureSet(all_merged, features = SARS.genes, col.name = "percent.viral", assay="Spatial")
# all_merged <- PercentageFeatureSet(all_merged, features = SARS.genes, col.name = "percent.viral", assay="RNA")

###### UMAP + Spatial DimPlots
plot2 <- SpatialFeaturePlot(all_merged, features="XIST", ncol = 1) + patchwork::plot_layout(ncol = 1)
####### FOR SOME REASON, INTEGRATION ADDS IMAGES ARBITRARILY
ggsave("integrated_SpatialDimPlot_XIST_TEST.pdf", plot = plot2, device = "pdf", width = 8, height = 45, units = "in")



##### Perform differential expression between seurat_clusters
dir.create("DE_seurat_clusters")
setwd("DE_seurat_clusters")
	## all_merged <- SCTransform(all_merged, method = "glmGamPoi", assay = "SCT", conserve.memory=TRUE, verbose = T)
# all_merged <- SCTransform(all_merged, method = "glmGamPoi", assay = "SCT", return.only.var.genes = FALSE, verbose = T)
DefaultAssay(all_merged) <- "SCT"
# all_merged <- NormalizeData(all_merged, normalization.method = "LogNormalize", scale.factor = 10000)
# all_merged <- ScaleData(all_merged, features = all.genes, assay="SCT")
# all_merged <- SCTransform(all_merged, method = "glmGamPoi", assay = "SCT", return.only.var.genes=FALSE, verbose = T)
var.genes <- all_merged@assays[["SCT"]]@var.features
# all_merged<- ScaleData(all_merged, features = var.genes, assay = "SCT", verbose = TRUE)

Idents(all_merged) <- "seurat_clusters"
#de_markers <- FindAllMarkers(all_merged, features = var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
#write.table(de_markers, "DEGs_by_seurat_clusters_pos-0.log2FC.txt", sep="\t")
de_markers <- FindAllMarkers(all_merged, features = var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_seurat_clusters_pos-0.693lnFC.txt", sep="\t")

cluster.averages <- AverageExpression(all_merged, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="SCT")
ggsave("SARS.genes_mc_heatmap_seurat_clusters.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="SCT")
ggsave("SARS.genes_FC_heatmap_seurat_clusters.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(all_merged, features = unique.top2)
png(file=paste0("seurat_clusters_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(all_merged, features = unique.top2)
png(file=paste0("seurat_clusters_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(all_merged, features = unique.top2)
png(file=paste0("seurat_clusters_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
## Slim down the all_mergeds to save space. Don't need to keep the scale.data
# all_merged <- DietSeurat(all_merged, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
setwd("..")



dir.create("immunology.clusters")
setwd("immunology.clusters")
Idents(all_merged) <- "seurat_clusters"
## use RNA for all DE analysis and plots
DefaultAssay(all_merged) <- "SCT"
# all_merged <- ScaleData(all_merged, features = all.genes, assay = "SCT")
all_merged <- SCTransform(all_merged, method = "glmGamPoi", assay = "SCT", return.only.var.genes = FALSE, verbose = T)

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
Idents(all_merged) <- "seurat_clusters"

# SARS.genes <- c("ORF1AB-1-1", "S", "ORF3A", "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF10")
p1 <- plot_density(all_merged, "CD4")
p2 <- FeaturePlot(all_merged, "CD8A")
p5 <- CombinePlots(plots = list(p1, p2))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("UMAP_KernelDensity_CD4-CD8A.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")

## 2 Density Plots
p3 <- plot_density(all_merged, c("CD4", "CD8A"))
p5 <- p3 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_CD4-CD8A_Joint.pdf", plot = p5, device = "pdf", width = 5, height = 7, units = "in")



## 2 Density Plots and a Joint density plot
p4 <- plot_density(all_merged, c("CD4", "CD8A"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_major.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(all_merged, c("TRGV9", "TRDV2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_gd.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(all_merged, c("TRAV10", "TRAJ18"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MAIT.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(all_merged, c("CD3G", "FCGR3A","NCAM1", "NCR1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_NK.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(all_merged, c("IFNG", "TBX21", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH1.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 16, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(all_merged, c("GATA3", "IL4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH2.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

## 2 Density Plots and a Joint density plot
p4 <- plot_density(all_merged, c("RORC", "IL17A", "IL17F", "IL21"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TH17.t_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 22, units = "in")




## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(all_merged, features = lymphoid.lineage)
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
feature.plot <- DotPlot(all_merged, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(all_merged, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

macrophage.lineage <- c("CD14", "CD68", "CSF1R", "ADGRE1", "SIGLEC1", "CD38", "GPR18", "FPR2", "MYC", "EGR2")
feature.plot <- DotPlot(all_merged, features = macrophage.lineage)
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

feature.plot <- DotPlot(all_merged, features = macrophage.lineage.markers.2)
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

#all_merged <- ScaleData(all_merged, features = all.genes, assay = "SCT")
cluster.averages <- AverageExpression(all_merged, assays="SCT", slot="counts", return.seurat = TRUE, features=t.lineages)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, slot="counts", assay="SCT")
ggsave("t.lineages_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, assay="SCT")
ggsave("t.lineages_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(all_merged, assays="SCT", slot="counts", return.seurat = TRUE, features=macrophage.lineage.markers.2)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, slot="counts", assay="SCT")
ggsave("macrophage.lineage.markers.2_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, assay="SCT")
ggsave("macrophage.lineage.markers.2_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(all_merged, assays="SCT", slot="counts", return.seurat = TRUE, features=lymphoid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, slot="counts", assay="SCT")
ggsave("lymphoid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, assay="SCT")
ggsave("lymphoid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(all_merged, assays="SCT", slot="counts", return.seurat = TRUE, features=myeloid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, slot="counts", assay="SCT")
ggsave("myeloid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, assay="SCT")
ggsave("myeloid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


############## Activation gating strategy https://www.nature.com/articles/s41590-021-01049-2/figures/10
## 2 Density Plots and a Joint density plot
p4 <- plot_density(all_merged, c("CD3D", "NCAM1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_TorNK_kernel-and-Joint.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(all_merged, c("CD4", "KLRB1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRB1-CD161.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(all_merged, c("CD4", "CD38"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD38.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
## 2 Density Plots and a Joint density plot
p4 <- plot_density(all_merged, c("CD4", "CD69"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-CD69.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("CD4", "KLRK1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv_kernel-and-KLRK1-NKG2D.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("CD4", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("FOXP3", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-FOXP3-CD4.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("KIT", "CD4"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Treg.cytokine_kernel-and-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("TNF", "TNFRSF8"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Tactiv.cytokine_kernel-and-TNF.TNFRSF8.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")

p4 <- plot_density(all_merged, c("CD14", "ITGAM"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("CD38", "FPR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_Macrophage_kernel-and-CD14-ITGAM.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("SIGLEC1", "FOLR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("SIGLEC1", "DDX3Y"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_HBC-fetalsex_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
	## Placenta-associated maternal macrophages PAMMS (FOLR2 neg, HLA positive); HBCs FOLR2+HLA- ; Thomas et al. Naomi McGovern, JEM 2020
p4 <- plot_density(all_merged, c("FOLR2", "HLA-DRA"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMMs_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM1_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("CD9", "CCR2"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_PAMM2_kernel.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("ITGAM", "CD14"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M0.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("CD80", "TNF"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("CD11B", "MYC"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("CD163", "MRC1"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2a.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("HLA-DRA", "IL6"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("CD14", "CD163"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_M2c.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("CD14", "CD1C"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p4 <- plot_density(all_merged, c("CD14", "THBD"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_KernelDensity_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
p1 <- plot_density(all_merged, "THBD")
p2 <- FeaturePlot(all_merged, "THBD")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(all_merged, "ITGB3")
p2 <- FeaturePlot(all_merged, "ITGB3")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_Megakaryocyte.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(all_merged, "CD1C")
p2 <- FeaturePlot(all_merged, "CD1C")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeloidDC1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(all_merged, "IL6")
p2 <- FeaturePlot(all_merged, "IL6")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2b.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(all_merged, "TNF")
p2 <- FeaturePlot(all_merged, "TNF")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M1.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(all_merged, "MYC")
p2 <- FeaturePlot(all_merged, "MYC")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_M2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(all_merged, "FOLR2")
p2 <- FeaturePlot(all_merged, "FOLR2")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_HBCs-FOLR2.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(all_merged, "HLA-DRA")
p2 <- FeaturePlot(all_merged, "HLA-DRA")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_PAMM-HLA-DRA.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p1 <- plot_density(all_merged, "CD33")
p2 <- FeaturePlot(all_merged, "CD33")
p5 <- CombinePlots(plots = list(p1, p2))
ggsave("UMAP_MyeliodProgenitor-M0.pdf", plot = p5, device = "pdf", width = 6, height = 3, units = "in")
p4 <- plot_density(all_merged, c("KIT", "CD33"), joint = TRUE)
p5 <- p4 + plot_layout(ncol = 1)
ggsave("UMAP_MyeliodProgenitor-M0-KIT.pdf", plot = p5, device = "pdf", width = 6, height = 11, units = "in")
# all_merged <- DietSeurat(all_merged, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
setwd("..")
rm(cluster.averages)
rm(library.averages.heatmap)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(p5)
rm(p6)
rm(p7)
rm(p8)
rm(plot2)
rm(umap.combined)
rm(feature.plot)
rm(de_markers)

setwd("..")
# all_merged <- DietSeurat(all_merged, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
# all_merged <- DietSeurat(all_merged, counts=TRUE, data=TRUE, scale.data = FALSE)





set.seed(seed=1)
setwd("/home/ebarrozo/visium/results/seurat_human_v2/integration-v3")
library(dplyr)
library(Seurat)	 # Seurat_4.0.3
library(sctransform)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(glmGamPoi)
library(ggpubr)
load("human-merged_v3.RData")


all_merged <- PercentageFeatureSet(all_merged, features = Spike.genes, col.name = "percent.Spike", assay="SCT")
all_merged <- PercentageFeatureSet(all_merged, features = HB.genes, col.name = "percent.Hemo", assay="SCT")
all_merged <- PercentageFeatureSet(all_merged, features = "DDX3Y", col.name = "percent.DDX3Y", assay="SCT")

SARS.receptors <- c("ACE2", "TMPRSS2", "CTSB", "CTSL", "HBA1", "HBB", "HBD", "HBZ", "HBQ1", "ICAM4", "DDX3Y", "SPP1")
SARS.genes.receptors.hb.fetalsex.mac <- union(SARS.genes, SARS.receptors)


Idents(all_merged) <- "orig.ident"
new.metadata <- c("M", "M", "F", "M", "F", "F", "F", "F", "F", "M", "M", "M", "M", "M", "F", "F", "M", "F", "F", "F", "M", "M", "M", "M", "M", "M", "M", "M", "M", "F", "F", "F", "M", "F", "M", "F", "F", "M", "M", "F", "M", "M", "M", "M", "M", "F", "F", "F", "F", "M", "F", "M", "M", "F", "F", "M", "F", "M", "M", "M", "M", "M", "M", "F", "F", "M", "M", "F", "F", "M", "F", "F", "F", "F", "F", "F", "F", "F", "F")
names(new.metadata) <- levels(all_merged)
all_merged <- RenameIdents(all_merged, new.metadata)
all_merged$FetalSex <- Idents(all_merged)


dir.create("rcpa.integration")
setwd("rcpa.integration")


Idents(all_merged) <- "etal"

all_merged <- SCTransform(all_merged, method = "glmGamPoi", assay = "SCT", return.only.var.genes=TRUE, do.scale=FALSE, conserve.memory = TRUE, verbose = TRUE)
all_merged <- FindVariableFeatures(all_merged)

all_merged.list <- SplitObject(all_merged, split.by = "etal")
reference.list <- all_merged.list[c("Tsang-2017", "Pique-Regi-2019", "Pique-Regi-2020", "Yang-2021", "Garcia-Flores-2022", "Barrozo-2022")]

reference.features <- SelectIntegrationFeatures(object.list = reference.list, nfeatures = 3000)
reference.list <- PrepSCTIntegration(object.list = reference.list, anchor.features = reference.features, verbose = T)

all_merged<- ScaleData(all_merged, features = reference.features, assay = "SCT", verbose = TRUE)
all_merged<- RunPCA(all_merged, features = reference.features , assay = "SCT", slot = "scale.data", verbose = TRUE)
all_merged.list <- SplitObject(all_merged, split.by = "etal")
rm("all_merged")
reference.list <- all_merged.list[c("Tsang-2017", "Pique-Regi-2019", "Pique-Regi-2020", "Yang-2021", "Garcia-Flores-2022", "Barrozo-2022")]
rm(all_merged.list)

reference.list <- PrepSCTIntegration(object.list = reference.list, anchor.features = reference.features, verbose = T)
all_merged.anchors <- FindIntegrationAnchors(object.list = reference.list, normalization.method = "SCT", reduction = "rpca", dims = 1:50, anchor.features = reference.features, verbose = T)
seurat.object <- IntegrateData(anchorset = all_merged.anchors, dims = 1:50)
rm(reference.list)
rm(all_merged.anchors)

# Use the 'integrated' assay for clustering and the SCT assay for differential expression
DefaultAssay(seurat.object) <- "integrated"
# Add all integrated genes as a list of features to call upon during clustering
integrated.genes <- rownames(seurat.object)
write.table(integrated.genes, "integrated.genes.txt", sep="\t")

seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE)
# save.image("human-rcpa-integrated_v3.RData")
# load("human-rcpa-integrated_v3.RData")
	# # save.image("human-rcpa-integrated_v1.RData")
############################################################# ############################################################# #############################################################
############################################################# Integrated: Dimension Reduction and Visualization #############################################################
############################################################# ############################################################# #############################################################

set.seed(seed=1)
setwd("/home/ebarrozo/visium/results/seurat_human_v2/integration-v3/rcpa.integration")
library(dplyr)
library(Seurat)	 # Seurat_4.0.3
library(sctransform)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(glmGamPoi)
library(ggpubr)
# load("human-rcpa-integrated_v1.RData")
rm(all_merged.anchors)

# seurat.object <- SCTransform(seurat.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes=TRUE, do.scale=FALSE, conserve.memory = TRUE, verbose = TRUE)
DefaultAssay(seurat.object) <- "SCT"
seurat.object <- PercentageFeatureSet(seurat.object, pattern = "^MT-", col.name = "percent.mt")
seurat.object <- PercentageFeatureSet(seurat.object, pattern = "^RP[SL]", col.name = "percent.ribo")
## Log normalize data and multiply by a factor of 10000
  ## required for cellcyclescoring
seurat.object<- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
# A list of human cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.gene
# Run cell cycle scoring 
seurat.object <- CellCycleScoring(seurat.object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
seurat.object$CC.Difference <- seurat.object$S.Score - seurat.object$G2M.Score

DefaultAssay(seurat.object) <- "integrated"
seurat.object <- ScaleData(seurat.object, features = integrated.genes, assay = "integrated", vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"))
		# Regressing out CC.Difference, percent.mt, percent.ribo
		#  |                                                                   |   0%
		# Error in qr.resid(qr = qr, y = data.expr[x, ]) : 
		#  'qr' and 'y' must have the same number of rows
seurat.object<- RunPCA(seurat.object, features = integrated.genes, assay = "integrated", slot = "scale.data")
seurat.object <- FindNeighbors(seurat.object, features = "integrated.genes", dims = 1:30)
seurat.object <- FindClusters(seurat.object, resolution = 0.6)
seurat.object <- RunUMAP(seurat.object, dims=1:30)
	### 41 clusters

seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

Idents(seurat.object) <- "seurat_clusters"
t1<-table(Idents(seurat.object))
write.table(t1, "seurat_clusters.counts.txt", sep="\t")
rm(t1)

# save.image("human-rcpa-integrated-umap_v3.RData")


set.seed(seed=1)
setwd("/home/ebarrozo/visium/results/seurat_human_v2/integration-v3/rcpa.integration")
library(dplyr)
library(Seurat)	 # Seurat_4.0.3
library(sctransform)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(glmGamPoi)
library(ggpubr)
load("human-rcpa-integrated-umap_v3.RData")



library(ggpubr)
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "etal", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='et al.')
p1
ggsave("UMAP_XL-etal.pdf", plot = p1, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Profile')
p2
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Site", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='Site')
ggsave("UMAP-Site.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Platform", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='Platform')
ggsave("UMAP-Platform.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Condition", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='Condition')
ggsave("UMAP-Condition.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")

umap.combined <- p1 + p2 + p3 + p4
ggsave("integrated_UMAP-FINAL.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("UMAP_XL.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")

p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "FetalSex", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='FetalSex')
ggsave("UMAP-FetalSex.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")
p4 <- FeaturePlot(seurat.object, features="DDX3Y", reduction = "umap", cols = mypal3, raster=FALSE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='DDX3Y')
ggsave("UMAP-DDX3Y.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")

#examine UMAPs with qc metrics 
p5 <- FeaturePlot(seurat.object, features = 'nFeature_RNA')
p6 <- FeaturePlot(seurat.object, features = 'nCount_RNA')
p7 <- FeaturePlot(seurat.object, features = 'percent.mt')
p8 <- FeaturePlot(seurat.object, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots-RNA.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")
#examine UMAPs with qc metrics 
p5 <- FeaturePlot(seurat.object, features = 'nFeature_Spatial'
p6 <- FeaturePlot(seurat.object, features = 'nCount_Spatial')
p7 <- FeaturePlot(seurat.object, features = 'percent.mt')
p8 <- FeaturePlot(seurat.object, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots-Spatial.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")
#examine UMAPs with qc metrics 
p5 <- DimPlot(seurat.object, reduction = "umap", group.by = "Phase", cols = mypal3, raster=FALSE) + labs(title = NULL, color='Cell Cycle Phase')
p6 <- FeaturePlot(seurat.object, features = 'percent.Hemo')
p7 <- FeaturePlot(seurat.object, features = 'percent.Spike')
p8 <- FeaturePlot(seurat.object, features = 'percent.viral')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots2.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")


Idents(seurat.object) <- "Platform"
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "Platform", label=TRUE)
ggsave("integrated_UMAP_splitby_platform-platform.pdf", plot = p2, device = "pdf", width = 9, height = 3, units = "in")
Idents(seurat.object) <- "Site"
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "Site", label=TRUE)
ggsave("integrated_UMAP_splitby_Site-Site.pdf", plot = p2, device = "pdf", width = 12, height = 3, units = "in")
Idents(seurat.object) <- "Condition"
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "Condition", label=TRUE)
ggsave("integrated_UMAP_splitby_Condition-Condition.pdf", plot = p2, device = "pdf", width = 12, height = 3, units = "in")
Idents(seurat.object) <- "Type"
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "etal", label=TRUE)
ggsave("integrated_UMAP_splitby_etal-Type.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
Idents(seurat.object) <- "Type"
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "Condition", label=TRUE)
ggsave("integrated_UMAP_splitby_Condition-Type.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
Idents(seurat.object) <- "Type"
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "Site", label=TRUE)
ggsave("integrated_UMAP_splitby_Site-Type.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")


###### UMAP + Spatial DimPlots
DefaultAssay(seurat.object) <- "SCT"
plot2 <- SpatialFeaturePlot(seurat.object, features="SIGLEC1", ncol = 1) + patchwork::plot_layout(ncol = 1)
####### FOR SOME REASON, INTEGRATION ADDS IMAGES ARBITRARILY
ggsave("integrated_SpatialDimPlot_SIGLEC1_TEST.pdf", plot = plot2, device = "pdf", width = 8, height = 45, units = "in")
	## determine the correct imageID by copying and pasting them here

DefaultAssay(seurat.object) <- "SCT"
Idents(seurat.object) <- "etal"
seurat.object <- ScaleData(seurat.object, features = SARS.genes.receptors.hb.fetalsex.mac, assay="SCT")
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes.receptors.hb.fetalsex.mac)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, slot="counts", assay="SCT")
ggsave("SARS.genes.receptors.hb.fetalsex.mac_mc_heatmap_etal.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, assay="SCT")
ggsave("SARS.genes.receptors.hb.fetalsex.mac_FC_heatmap_etal.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

Idents(seurat.object) <- "Type"
seurat.object <- ScaleData(seurat.object, features = SARS.genes.receptors.hb.fetalsex.mac, assay="SCT")
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes.receptors.hb.fetalsex.mac)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, slot="counts", assay="SCT")
ggsave("SARS.genes.receptors.hb.fetalsex.mac_mc_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, assay="SCT")
ggsave("SARS.genes.receptors.hb.fetalsex.mac_FC_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

Idents(seurat.object) <- "seurat_clusters"
seurat.object <- ScaleData(seurat.object, features = SARS.genes.receptors.hb.fetalsex.mac, assay="SCT")
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes.receptors.hb.fetalsex.mac)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, slot="counts", assay="SCT")
ggsave("SARS.genes.receptors.hb.fetalsex.mac_mc_heatmap_seurat_clusters.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, assay="SCT")
ggsave("SARS.genes.receptors.hb.fetalsex.mac_FC_heatmap_seurat_clusters", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


dir.create("immunology.clusters")
setwd("immunology.clusters")
Idents(seurat.object) <- "seurat_clusters"
## use RNA for all DE analysis and plots
DefaultAssay(seurat.object) <- "SCT"
# seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "SCT")
seurat.object <- ScaleData(seurat.object, features = all.genes, assay="SCT")

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

#seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "SCT")
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=t.lineages)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, slot="counts", assay="SCT")
ggsave("t.lineages_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = t.lineages, raster = FALSE, assay="SCT")
ggsave("t.lineages_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=macrophage.lineage.markers.2)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, slot="counts", assay="SCT")
ggsave("macrophage.lineage.markers.2_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = macrophage.lineage.markers.2, raster = FALSE, assay="SCT")
ggsave("macrophage.lineage.markers.2_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=lymphoid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, slot="counts", assay="SCT")
ggsave("lymphoid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = lymphoid.lineage, raster = FALSE, assay="SCT")
ggsave("lymphoid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=myeloid.lineage)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, slot="counts", assay="SCT")
ggsave("myeloid.lineage_mc_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = myeloid.lineage, raster = FALSE, assay="SCT")
ggsave("myeloid.lineage_FC_heatmap_library.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


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
# seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
setwd("..")

rm(cluster.averages)
rm(library.averages.heatmap)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(p5)
rm(p6)
rm(p7)
rm(p8)
rm(plot2)
rm(umap.combined)
rm(feature.plot)
rm(de_markers)

##### Perform differential expression between seurat_clusters
dir.create("DE_seurat_clusters")
setwd("DE_seurat_clusters")
	## all_merged <- SCTransform(all_merged, method = "glmGamPoi", assay = "SCT", conserve.memory=TRUE, verbose = T)
# seurat.object <- SCTransform(seurat.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes = FALSE, verbose = T)
DefaultAssay(seurat.object) <- "SCT"
# seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
# seurat.object <- ScaleData(seurat.object, features = all.genes, assay="SCT")
seurat.object <- SCTransform(seurat.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes=TRUE, conserve.memory=TRUE, verbose = T)

var.genes <- seurat.object@assays[["SCT"]]@var.features
combined.var.genes <- union(var.genes, SARS.genes.receptors.hb.fetalsex.mac)
seurat.object <- NormalizeData(seurat.object, features = combined.var.genes, assay="SCT")
seurat.object <- ScaleData(seurat.object, features = combined.var.genes, assay="SCT")
Idents(seurat.object) <- "seurat_clusters"
#de_markers <- FindAllMarkers(seurat.object, features = var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
#write.table(de_markers, "DEGs_by_seurat_clusters_pos-0.log2FC.txt", sep="\t")
# de_markers <- FindAllMarkers(seurat.object, features = combined.var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
de_markers <- FindAllMarkers(seurat.object, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_seurat_clusters_pos-0.693lnFC.txt", sep="\t")

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="SCT")
ggsave("SARS.genes_mc_heatmap_seurat_clusters.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="SCT")
ggsave("SARS.genes_FC_heatmap_seurat_clusters.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat_clusters_top2.marker-DotPlot.png"),res=300, width=7500, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="SCT")
ggsave("top2_mc_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="SCT")
ggsave("top2_FC_heatmap.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat_clusters_top3.marker-DotPlot.png"),res=300, width=9000, height=3000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat_clusters_top5.marker-DotPlot.png"),res=300, width=10000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()


## Slim down the seurat.objects to save space. Don't need to keep the scale.data
# seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
setwd("..")

##### Perform differential expression between Type
dir.create("DE_Type")
setwd("DE_Type")
#DefaultAssay(seurat.object) <- "SCT"
#seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
#seurat.object <- ScaleData(seurat.object, features = all.genes, assay="SCT")
Idents(seurat.object) <- "Type"
# de_markers <- FindAllMarkers(seurat.object, features = combined.var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
de_markers <- FindAllMarkers(seurat.object, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)

write.table(de_markers, "DEGs_by_Type_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="SCT")
ggsave("SARS.genes_mc_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="SCT")
ggsave("SARS.genes_FC_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Type_top2.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Type_top3.marker-DotPlot.png"),res=300, width=4500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Type_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = p_val_adj)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Type_top5.marker-DotPlot-p_val_adj.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="SCT")
ggsave("top5_mc_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="SCT")
ggsave("top5_FC_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


setwd("..")


##### Perform differential expression between Site
dir.create("DE_Site")
setwd("DE_Site")
DefaultAssay(seurat.object) <- "SCT"
# seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
# seurat.object <- ScaleData(seurat.object, features = all.genes, assay="SCT")
Idents(seurat.object) <- "Site"
# de_markers <- FindAllMarkers(seurat.object, features = combined.var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
de_markers <- FindAllMarkers(seurat.object, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)

write.table(de_markers, "DEGs_by_Site_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes)
####Use this to determine the sum of mean counts for macrophage.lineage.markers.2 transcripts per cluster and rank macrophage.lineage.markers.2 transcripts based on overall macrophage.lineage.markers.2 transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="SCT")
ggsave("SARS.genes_mc_heatmap_Site.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for macrophage.lineage.markers.2 transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="SCT")
ggsave("SARS.genes_FC_heatmap_Site.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Site_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Site_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Site_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="SCT")
ggsave("top5_mc_heatmap_Condition.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="SCT")
ggsave("top5_FC_heatmap_Condition.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


setwd("..")

## Slim down the seurat.objects to save space. Don't need to keep the scale.data
# seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")


##### Perform differential expression between Condition
dir.create("DE_Condition")
setwd("DE_Condition")
#DefaultAssay(seurat.object) <- "SCT"
#seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
#seurat.object <- ScaleData(seurat.object, features = all.genes, assay="SCT")
Idents(seurat.object) <- "Condition"
# de_markers <- FindAllMarkers(seurat.object, features = combined.var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
de_markers <- FindAllMarkers(seurat.object, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)

write.table(de_markers, "DEGs_by_Condition_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="SCT")
ggsave("SARS.genes_mc_heatmap_Condition.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="SCT")
ggsave("SARS.genes_FC_heatmap_Condition.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Condition_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Condition_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Condition_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = p_val_adj)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("Condition_top5.marker-DotPlot-p_val_adj.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="SCT")
ggsave("top5_mc_heatmap_Condition.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="SCT")
ggsave("top5_FC_heatmap_Condition.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

Idents(seurat.object) <- "Type"
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="SCT")
ggsave("top5_mc_heatmap_Condition-Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="SCT")
ggsave("top5_FC_heatmap_Condition-Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")



setwd("..")

##### Perform differential expression between FetalSex
dir.create("DE_FetalSex")
setwd("DE_FetalSex")
#DefaultAssay(seurat.object) <- "SCT"
#seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
#seurat.object <- ScaleData(seurat.object, features = all.genes, assay="SCT")
Idents(seurat.object) <- "FetalSex"
# de_markers <- FindAllMarkers(seurat.object, features = combined.var.genes, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
de_markers <- FindAllMarkers(seurat.object, assay = "SCT", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "DEGs_by_FetalSex_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, slot="counts", assay="SCT")
ggsave("SARS.genes_mc_heatmap_FetalSex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes, raster = FALSE, assay="SCT")
ggsave("SARS.genes_FC_heatmap_FetalSex.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("FetalSex_top2.marker-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("FetalSex_top3.marker-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("FetalSex_top5.marker-DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = p_val_adj)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("FetalSex_top5.marker-DotPlot-p_val_adj.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, slot="counts", assay="SCT")
ggsave("top5_mc_heatmap_Condition.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = FALSE, assay="SCT")
ggsave("top5_FC_heatmap_Condition.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


setwd("..")

seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

dir.create("gene.set.analysis")
setwd("gene.set.analysis")

dir.create("hmo.type")
setwd("hmo.type")
top3k <- seurat.object@assays[["SCT"]]@var.features
DefaultAssay(seurat.object) <- "SCT"
# all.genes <- rownames(seurat.object@assays[["SCT"]])
seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "SCT")

Idents(seurat.object) <- "Type"

## HMO secreter status by FUT2 and Lewis by FUT3 but not HMO specific https://www.frontiersin.org/articles/10.3389/fnut.2020.574459/full. Also LALBA from Max's converstation with Lars.
hmo.markers <- c("FUT2", "FUT3", "LALBA", "DDX3Y", "XIST")

feature.plot <- FeaturePlot(seurat.object, features = hmo.markers, pt.size = 0.0001)
ggsave("hmo.markers-FeaturePlots.pdf", plot = feature.plot, device = "pdf", width = 8, height = 15, units = "in")
biomarkers.heatmap <- DoHeatmap(seurat.object, features = hmo.markers, raster = TRUE)
ggsave("hmo.markers-logfc_heatmap.pdf", plot = biomarkers.heatmap, device = "pdf", width = 17, height = 5, units = "in")
biomarkers.heatmap <- DoHeatmap(seurat.object, features = hmo.markers, slot="counts", raster = TRUE)
ggsave("hmo.markers-counts_heatmap.pdf", plot = biomarkers.heatmap, device = "pdf", width = 17, height = 5, units = "in")
feature.plot <- DotPlot(seurat.object, features = hmo.markers)
png(file=paste0("hmo.markers-DotPlot.png"),
                res=300, 
                width=2500, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="HMO Markers") + scale_y_discrete(name ="Profile")
dev.off()

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=hmo.markers)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = hmo.markers, raster = FALSE, assay="SCT")
ggsave("hmo.markers_FC_heatmap_library.pdf", plot = library.averages.heatmap, device = "pdf", width = 12, height = 14, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = hmo.markers, raster = FALSE, slot="counts", assay="SCT")
ggsave("hmo.markers_mc_heatmap_library.pdf", plot = library.averages.heatmap, device = "pdf", width = 12, height = 14, units = "in")


######## HMO enzymes
## From https://www.biorxiv.org/content/10.1101/2020.09.02.278663v1.full, Supp. Table 1, filtered for Candidate for Further Consideration: yes and not determined
hmo.enzymes <- c("B3GNT2", "B3GNT3", "B3GNT6", "B3GNT8", "B3GNT9", "FUT2", "FUT3", "FUT4", "FUT5", "FUT6", "FUT9", "FUT10", "FUT11", "ST3GAL1", "ST3GAL3", "ST3GAL5", "ST6GAL1", "B3GALT4", "B3GALT5", "B4GALT1", "B4GALT2", "B4GALT3", "B4GALT4", "B4GALT5", "GCNT1", "GCNT2A", "GCNT3", "FUT3", "FUT5", "ST6GALNAC2", "ST6GALNAC3", "ST6GALNAC4", "ST6GALNAC6", "DDX3Y", "XIST")
hmo.enzymes <- unique(hmo.enzymes)
hmo.top.markers <- c(top2, hmo.enzymes)
hmo.top.markers <- unique(hmo.top.markers)

#feature.plot <- FeaturePlot(seurat.object, features = hmo.enzymes, pt.size = 0.0001)
#ggsave("hmo.enzymes-FeaturePlots.pdf", plot = feature.plot, device = "pdf", width = 8, height = 15, units = "in")
biomarkers.heatmap <- DoHeatmap(seurat.object, features = hmo.enzymes, raster = FALSE)
ggsave("hmo.enzymes-logfc_heatmap.pdf", plot = biomarkers.heatmap, device = "pdf", width = 20, height = 5, units = "in")
biomarkers.heatmap <- DoHeatmap(seurat.object, features = hmo.enzymes, slot="counts", raster = FALSE)
ggsave("hmo.enzymes-counts_heatmap.pdf", plot = biomarkers.heatmap, device = "pdf", width = 20, height = 5, units = "in")

feature.plot <- DotPlot(seurat.object, features = hmo.enzymes)
png(file=paste0("hmo.enzymes-DotPlot.png"),
                res=300, 
                width=2500, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="HMO Enzymes") + scale_y_discrete(name ="Profile")
dev.off()

cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=hmo.enzymes)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = hmo.enzymes, raster = FALSE, assay="SCT", label=FALSE)
ggsave("hmo.enzymes_FC_heatmap_library.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = hmo.enzymes, raster = FALSE, slot="counts", assay="SCT", label=FALSE)
ggsave("hmo.enzymes_mc_heatmap_library.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

setwd("..")
setwd("..")











###### UMAP + Spatial DimPlots
# plot2 <- SpatialFeaturePlot(seurat.object, features="XIST", ncol = 1) + patchwork::plot_layout(ncol = 1)
####### FOR SOME REASON, INTEGRATION ADDS IMAGES ARBITRARILY
# ggsave("integrated_SpatialDimPlot_XIST_TEST.pdf", plot = plot2, device = "pdf", width = 8, height = 45, units = "in")
	## determine the correct imageID by copying and pasting them here
## copy and paste the title of image that works S01.1.15 	S03.1.13.13 	S04.2.14.14		S04.3.15

##############################################################################################
  # Metadata
      ## uninfected: S01, S02, S03, S24, S25, S26
      ## SARS-CoV-2 asymptomatic: S15, S17, S18, S19
      ## SARS-CoV-2 symptomatic: S16, S20, S21, S22, S23a, S23b
        ## SARS-CoV-2 perivillous fibrinoids: S23a and S23b
        ## Noted areas of inflammation: S18
      ## Samples S01, S02, and S03 were from the sample placenta from distinct regions including chorionic villi, decidua, and a chorioamiotic membrane roll
      ## Samples 15-26 were sampled from the placental parenchyme with various regions of villi and some decidua
      ## Samples 23a and 23b were from the same placenta punch, selected for two sections due to the presence of perivillious fibrinoid pathology identified in QC H&E slide

S01.21	S01.1.20	S03.1.5	S04.2.5		S15.3.5	S16.4.5	S17.4.5 S18.4.5	 S19.4.5 S20.4.5


Idents(seurat.object) <- "Type"
###### UMAP + Spatial DimPlots --- labeled by Type      ----- large dot sizes	-- color schemed
p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE, cols = mypal3)
p1
image.list <- c("slice1.5")
p2 <- SpatialDimPlot(seurat.object, images=image.list, pt.size.factor = 3, cols = mypal3)
p2
image.list <- c("slice1.1")
p3 <- SpatialDimPlot(seurat.object, images=image.list, cols = mypal3, pt.size.factor = 3)
p3
image.list <- c("slice1.2.1")
p4 <- SpatialDimPlot(seurat.object, images=image.list, cols = mypal3, pt.size.factor = 3)
p4
image.list <- c("slice1.3.3")
p5 <- SpatialDimPlot(seurat.object, images=image.list, cols = mypal3, pt.size.factor = 3)
p5
	######### Change sliceID for each group
p6 <- wrap_plots(p1, p2, p3, p4, p5, ncol=7)
ggsave("Type_annotated_UMAP_SpatialDimPlots-cropped-Final-Uninfected.pdf", plot = p6, device = "pdf", width = 28, height = 4, units = "in")







p6 <- wrap_plots(p1, p2, p3, p4, p5, ncol=5)
ggsave("Type_annotated_UMAP_SpatialDimPlots-cropped-Final-SARS-CoV-2.ax.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")







p6 <- wrap_plots(p1, p2, p3, p4, p5, ncol=7)
ggsave("Type_annotated_UMAP_SpatialDimPlots-cropped-Final-SARS-CoV-2.sx.pdf", plot = p6, device = "pdf", width = 28, height = 4, units = "in")




## SpatialPlots with all samples showing a particular gene
Idents(seurat.object) <- "SampleID"
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=2, cols = mypal3, pt.size.factor = 4)
p5
ggsave("C3_SpatialPlots-cropped-Final.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=2, alpha = c(0.1, 6), cols = mypal3, pt.size.factor = 4)
p5
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=2, max.cutoff=6, cols = mypal3, pt.size.factor = 4)
p5
ggsave("C3_SpatialPlots-cropped-same.scale-Final.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

image.list <- c("slice1.5")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, alpha = c(0.1, 1), max.cutoff = "q95", cols = mypal3, pt.size.factor = 4)
p5

image.list <- c("slice1.1")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, pt.size.factor = 4, max.cutoff = "7") + ggplot2::scale_fill_binned(limits = c(0.0,7), breaks = c(0.0, 3.5, 7.0))
p5

image.list <- c("slice1.1","slice1.2.1","slice1.5", "slice1.3.3")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=2, pt.size.factor = 4, combine=T, max.cutoff=2) 
ggsave("C3_SpatialPlots-cropped-Final.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

##############################################################################################





p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("S01.1.15")
p2 <- SpatialDimPlot(seurat.object, images=image.list)
p2
image.list <- c("S03.1.13.13")
p3 <- SpatialDimPlot(seurat.object, images=image.list)
p3
image.list <- c("S04.2.14.14")
p4 <- SpatialDimPlot(seurat.object, images=image.list)
p4
image.list <- c("S15.3.15.15")
p5 <- SpatialDimPlot(seurat.object, images=image.list)
p5
	## Check to confirm they are in the correct order
p6 <- wrap_plots(p1, p2, p3, p4, p5, ncol=5)
ggsave("integrated_UMAP_SpatialDimPlots-cropped.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")

p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("S01.1.15")
p2 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p2
image.list <- c("S03.1.13")
p3 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p3
image.list <- c("S04.2.14")
p4 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p4
image.list <- c("S15.3.15")
p5 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p5
p6 <- wrap_plots(p1, p2, p3, p4, p5, ncol=5)
ggsave("integrated_UMAP_SpatialDimPlots-uncropped.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")

Idents(seurat.object) <- "Type"

p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("S01.1.15")
p2 <- SpatialDimPlot(seurat.object, images=image.list)
p2
image.list <- c("S03.1.13")
p3 <- SpatialDimPlot(seurat.object, images=image.list)
p3
image.list <- c("S04.2.14")
p4 <- SpatialDimPlot(seurat.object, images=image.list)
p4
image.list <- c("S15.3.15")
p5 <- SpatialDimPlot(seurat.object, images=image.list)
p5
    ## Check to confirm they are in the correct order
p6 <- wrap_plots(p1, p2, p3, p4, p5, ncol=5)
ggsave("integrated_UMAP_SpatialDimPlots-cropped_type.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")

p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
p1
image.list <- c("S01.1.15")
p2 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p2
image.list <- c("S03.1.13")
p3 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p3
image.list <- c("S04.2.14")
p4 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p4
image.list <- c("S15.3.15")
p5 <- SpatialDimPlot(seurat.object, images=image.list, crop = FALSE)
p5
p6 <- wrap_plots(p1, p2, p3, p4, p5, ncol=5)
ggsave("integrated_UMAP_SpatialDimPlots-uncropped_type.pdf", plot = p6, device = "pdf", width = 20, height = 4, units = "in")


Idents(seurat.object) <- "seurat_clusters"


## Spatial DimPlots split.by clusters
image.list <- c("S04.2.14")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S07.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S01.1.15")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S09.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S03.1.13")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S12.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S15.3.15")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S13.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


image.list <- c("S04.2.14")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(10, 7, 17)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S07_macs.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S01.1.15")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(10, 7, 17)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S09_macs.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S03.1.13")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(10, 7, 17)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S12_macs.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S15.3.15")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(10, 7, 17)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S13_macs.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")


image.list <- c("S04.2.14")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0,1,3,4,5,6,8,9,11,13,14)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S07_NK.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S01.1.15")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0,1,3,4,5,6,8,9,11,13,14)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S09_NK.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S03.1.13")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0,1,3,4,5,6,8,9,11,13,14)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S12_NK.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S15.3.15")
p3<- SpatialDimPlot(seurat.object, images=image.list, cells.highlight = CellsByIdentities(object = seurat.object, idents = c(0,1,3,4,5,6,8,9,11,13,14)), facet.highlight = TRUE, ncol = 5, pt.size.factor = 4.0)
ggsave("integrated_UMAP_Spatial_SpatialDimPlots_split.by.cluster_S13_NK.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")




rm(umap.combined)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(p5)
rm(p6)
rm(p7)
rm(p8)
rm(p9)
rm(p10)
rm(image.list)

seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
# save.image("human_spatial_data-tsang-saha-integrated-annotated_v1.RData")










############################################################# ############################################################# #############################################################
############################################################# subset infected datasets only #############################################################
############################################################# ############################################################# #############################################################

setwd("/home/ebarrozo/visium/results/seurat_human_v2/integration-v3")
# load("human-integrated-umap_v3.RData")
dir.create("subset.SARSCOV2")
setwd("subset.SARSCOV2")

Idents(seurat.object) <- "Condition"
levels(seurat.object)
ncol(seurat.object)
	# 45413

seurat.object2 <- subset(x = seurat.object, idents = c("SARS-CoV-2.ax", "SARS-CoV-2.sx", "SARS-CoV-2"))
ncol(seurat.object2)
	# 20419;; 16511 tsang et al 2020
############################################################# ############################################################# #############################################################
############################################################# Integrated: Dimension Reduction and Visualization #############################################################
############################################################# ############################################################# #############################################################

seurat.object2 <- ScaleData(seurat.object2, features = integrated.genes, assay = "integrated", vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"))
seurat.object2<- RunPCA(seurat.object2, features = integrated.genes, assay = "integrated", slot = "scale.data")

seurat.object2 <- FindNeighbors(seurat.object2, features = "integrated.genes", dims = 1:30)
seurat.object2 <- FindClusters(seurat.object2, resolution = 0.6)
seurat.object2 <- RunUMAP(seurat.object2, dims=1:30)
	## clusters 0 thru 16
library(ggpubr)
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)


p1 <- DimPlot(seurat.object2, reduction = "umap", group.by = "etal", cols = mypal3, raster=FALSE) + labs(title = NULL, color='et al.')
ggsave("UMAP_XL-etal.pdf", plot = p1, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(seurat.object2, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Profile')
p2
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p3 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Site", cols = mypal3, raster=FALSE) + labs(title = NULL, color='Site')
ggsave("UMAP-Site.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p3 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Platform", cols = mypal3, raster=FALSE) + labs(title = NULL, color='Platform')
ggsave("UMAP-Platform.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p4 <- DimPlot(seurat.object2, reduction = "umap", group.by = "Condition", cols = mypal3, raster=FALSE) + labs(title = NULL, color='Condition')
ggsave("UMAP-Condition.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")

umap.combined <- p1 + p2 + p3 + p4
ggsave("integrated_UMAP-FINAL.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("UMAP_XL.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")








Idents(seurat.object) <- "orig.ident"
seurat.object <- subset(seurat.object, idents=("S01", "S03", "S04", "S15"))
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
p9 <- FeaturePlot(seurat.object, features = 'nFeature_RNA')
p10 <- FeaturePlot(seurat.object, features = 'nCount_RNA')
umap.combined <- CombinePlots(plots = list(p5, p9, p6, p10, p7, p8), ncol=2)
ggsave("integrated_UMAP_QCmetricsFeaturePlots.pdf", plot = umap.combined, device = "pdf", width = 8, height = 10, units = "in")


###### UMAP + Spatial DimPlots
## copy and paste the title of image that works S01.1.15 	S03.1.13.13 	S04.2.14.14		S04.3.15


rm(umap.combined)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(p5)
rm(p6)
rm(p7)
rm(p8)
rm(p9)
rm(p10)
rm(image.list)

setwd("..")



############################################################# ############################################################# #############################################################
############################################################# VISUALISATION OF TOP MARKERS #############################################################
############################################################# ############################################################# #############################################################


dir.create("Top.Markers")
setwd("Top.Markers")
## SpatialPlots with all samples showing a particular gene
image.list <- c("S04.2.14","S01.1.15","S03.1.13", "S15.3.15")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=2)
p5
ggsave("C3_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("S04.2.14","S01.1.15","S03.1.13", "S15.3.15")
p5 <- SpatialPlot(seurat.object, features = "C3", images=image.list, ncol=2, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("C3_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("S04.2.14","S01.1.15","S03.1.13", "S15.3.15")
p5 <- SpatialPlot(seurat.object, features = "KISS1", images=image.list, ncol=2)
p5
ggsave("KISS1_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("S04.2.14","S01.1.15","S03.1.13", "S15.3.15")
p5 <- SpatialPlot(seurat.object, features = "KISS1", images=image.list, ncol=2, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("KISS1_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("S04.2.14","S01.1.15","S03.1.13", "S15.3.15")
p5 <- SpatialPlot(seurat.object, features = "GDF15", images=image.list, ncol=2)
p5
ggsave("GDF15_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("S04.2.14","S01.1.15","S03.1.13", "S15.3.15")
p5 <- SpatialPlot(seurat.object, features = "GDF15", images=image.list, ncol=2, min.cutoff = 0, max.cutoff = 6)
p5
ggsave("GDF15_SpatialPlots-cropped-same.scale.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


## SpatialPlots with all samples showing a particular genes # DAF is CD55, inhibits complement
image.list <- c("S04.2.14","S01.1.15","S03.1.13", "S15.3.15")
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
image.list <- c("S04.2.14","S01.1.15","S03.1.13", "S15.3.15")
mock.markers <- c("LUM", "DKK1", "IGFBP1", "FN1", "PRG2")
p5 <- SpatialPlot(seurat.object, features = mock.markers, images=image.list, ncol=2)
p5
ggsave("Mock.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")

## SpatialPlots with all samples showing a particular gene
image.list <- c("S04.2.14","S01.1.15","S03.1.13", "S15.3.15")
marker.list <- c("PSG7", "CSH2", "KISS1", "PSG1", "GDF15")
p5 <- SpatialPlot(seurat.object, features = marker.list, images=image.list, ncol=2)
p5
ggsave("SARS-CoV-2.markers_SpatialPlots-cropped.pdf", plot = p5, device = "pdf", width = 10, height = 16, units = "in")


setwd("..")

############################################################# ############################################################# #############################################################
############################################################# DE Analysis #############################################################
############################################################# ############################################################# #############################################################
		## Unresolved error using FindAllMarkers. None of the assays are working. 
		## Unresolved error using FindAllMarkers. None of the assays are working. 
		## Unresolved error using FindAllMarkers. None of the assays are working. 
> de_markers <- FindAllMarkers(seurat.object, features = var.combined, assay = "RNA", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
Calculating cluster 0
Calculating cluster 1
Calculating cluster 2
Calculating cluster 3
Calculating cluster 4
Calculating cluster 5
Calculating cluster 6
Calculating cluster 7
Calculating cluster 8
Calculating cluster 9
Calculating cluster 10
Calculating cluster 11
Calculating cluster 12
Calculating cluster 13
Calculating cluster 14
Calculating cluster 15
Warning: No DE genes identified
Warning: The following tests were not performed: 
Warning: When testing 0 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 1 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 2 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 3 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 4 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 5 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 6 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 7 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 8 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 9 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 10 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 11 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 12 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 13 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 14 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
Warning: When testing 15 versus all:
	error in evaluating the argument 'x' in selecting a method for function 'rowSums': invalid character indexing
> 
		## Unresolved error using FindAllMarkers. None of the assays are working. 
		## Unresolved error using FindAllMarkers. None of the assays are working. 
		## Unresolved error using FindAllMarkers. None of the assays are working. 


seurat.object <- SCTransform(seurat.object, method = "glmGamPoi", assay="SCT", new.assay.name = "DE", do.scale=TRUE, return.only.var.genes= F, verbose = T)

DefaultAssay(seurat.object) <- "DE"
DE.var.combined <- seurat.object@assays[["DE"]]@var.features
##### Perform differential expression between clusters using the raw data
DefaultAssay(seurat.object) <- "Spatial"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = var.combined, assay="Spatial")

DefaultAssay(seurat.object) <- "RNA"
seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.object <- ScaleData(seurat.object, features = var.combined, assay="RNA")

Idents(seurat.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(seurat.object, features = var.combined, assay = "DE", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.693)
write.table(de_markers, "Spatial_DEGs_byclusters_pos-0.693lnFC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
View(top5)
	## make sure there are some transcripts that passed the lnFC(>0.693)
		## all clusters have sig. genes
			# de_markers <- FindAllMarkers(seurat.object, features = var.combined, assay = "Spatial", only.pos = TRUE, min.pct = 0.10, logfc.threshold =  0.2)
			# top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
			# View(top5)
			# write.table(de_markers, "integrated_DEGs_byclusters_pos-0.log2FC.txt", sep="\t")
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

image.list <- c("slice1.5")
p2 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p2
ggsave("top.markers_SpatialFeaturePlots-cropped-S07.pdf", plot = p2, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("S01.1.15")
p3 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p3
ggsave("top.markers_SpatialFeaturePlots-cropped-S09.pdf", plot = p3, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.2.1")
p4 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=7)
p4
ggsave("top.markers_SpatialFeaturePlots-cropped-S12.pdf", plot = p4, device = "pdf", width = 8, height = 10, units = "in")
image.list <- c("slice1.3.3")
p5 <- SpatialFeaturePlot(seurat.object, features = unique.top2, images=image.list, ncol=7)
p5
ggsave("top.markers_SpatialFeaturePlots-cropped-S13.pdf", plot = p5, device = "pdf", width = 8, height = 10, units = "in")


top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)

image.list <- c("slice1.5")
p2 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p2
ggsave("top.markers_SpatialPlots-cropped-S07.pdf", plot = p2, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("S01.1.15")
p3 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=5)
p3
ggsave("top.markers_SpatialPlots-cropped-S09.pdf", plot = p3, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.2.1")
p4 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=7)
p4
ggsave("top.markers_SpatialPlots-cropped-S12.pdf", plot = p4, device = "pdf", width = 10, height = 8, units = "in")
image.list <- c("slice1.3.3")
p5 <- SpatialPlot(seurat.object, features = unique.top2, images=image.list, ncol=7)
p5
ggsave("top.markers_SpatialPlots-cropped-S13.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")



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
############################################################# Gene set analysis #############################################################
############################################################# ############################################################# #############################################################
Idents(seurat.object) <- "Type"
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

feature.plot <- DotPlot(seurat.object, features = ZIKV.genes)
png(file=paste0("ZIKV.genes-DotPlot.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="ZIKV Transcripts") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
library.averages.heatmap <- DoHeatmap(seurat.object, features = ZIKV.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("ZIKV.genes_seurat.clusters_heatmap_counts.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = ZIKV.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("ZIKV.genes_seurat.clusters_heatmap_logfc.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
rm(library.averages.heatmap)

Idents(seurat.object) <- "Type"
seurat.object <- ScaleData(seurat.object, features = all.genes, assay = "Spatial")
## canonical markers for lymphoid cells https://docs.abcam.com/pdf/immunology/lymphoid_cell_subtype.pdf
#lymphoid.lineage <- c("CD34", "CD117", "CD271", "CD2", "CD3", "CD8A", "CD8B", "CD4", "CD25", "CD56", "CD30", "CD19", "CD20", "CD138", "CD304")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
lymphoid.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1")

feature.plot <- DotPlot(seurat.object, features = lymphoid.lineage)
png(file=paste0("lymphoid.lineage-DotPlot_subtype.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Lymphoid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

## canonical markers for myeloid cells https://docs.abcam.com/pdf/immunology/myeloid_cell_subtype.pdf
# myeloid.lineage <- c("CD34", "CD117", "CD271", "CD33", "CD71", "CD61", "CD123", "CD44", "CD15", "CD16", "CD11b", "CD14", "CD1c", "CD141")
    ## find all.genes.txt and convert common names to gene names listed in features, confirm with genecard.com
myeloid.lineage <- c("CD34", "KIT", "NGFR", "CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
feature.plot <- DotPlot(seurat.object, features = myeloid.lineage)
png(file=paste0("myeloid.lineage-DotPlot_subtype.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Myeloid Lineage Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

immunophenotype.lineage <- c("CD34", "KIT", "NGFR", "CD2", "CD3D","CD3E","CD3G", "CD8A", "CD8B", "CD8B2","CD4", "IL2RA", "NCAM1", "TNFRSF8", "CD19", "MS4A1", "SDC1", "NRP1","CD33", "TFRC", "ITGB3", "IL3RA", "CD44", "FUT4", "FCGR3A", "ITGAM", "CD14", "CD1C", "THBD")
immunophenotype.lineage <- unique(immunophenotype.lineage)
feature.plot <- DotPlot(seurat.object, features = immunophenotype.lineage)
png(file=paste0("immunophenotype.lineage-DotPlot_subtype.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Immunephenotype Markers") + scale_y_discrete(name ="UMAP Cluster")
dev.off()

feature.plot <- DotPlot(seurat.object, features = ZIKV.genes)
png(file=paste0("ZIKV.genes-DotPlot_subtype.png"),
                res=300, 
                width=3000, 
                height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="ZIKV Transcripts") + scale_y_discrete(name ="UMAP Cluster")
dev.off()
library.averages.heatmap <- DoHeatmap(seurat.object, features = ZIKV.genes, raster = FALSE, slot="counts", assay="Spatial", label=FALSE)
ggsave("ZIKV.genes_seurat.clusters_heatmap_counts_subtype.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(seurat.object, features = ZIKV.genes, raster = FALSE, slot="scale.data", assay="Spatial", label=FALSE)
ggsave("ZIKV.genes_seurat.clusters_heatmap_logfc_subtype.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
rm(library.averages.heatmap)

############################################################# ############################################################# #############################################################
seurat.object2 <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")
seurat.object <- seurat.object2
rm(seurat.object2)


#Save tables with cell counts per library and per cluster
dir.create("cell.counts")
setwd("cell.counts")
# Save a table with cell counts per library
Idents(seurat.object) <- "orig.ident"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "sample.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(seurat.object) <- "seurat_clusters"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "cluster.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(seurat.object) <- "Subtype"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "Subtype.cellcounts.txt", sep="\t")

# Save a table with cell counts per cluster
Idents(seurat.object) <- "Type"
unfiltered.count <- table(Idents(seurat.object))
write.table(unfiltered.count, "Type.cellcounts.txt", sep="\t")


dir.create("Type.cluster")
setwd("Type.cluster")

# "Trophoblast", "Fibroblast", "Stromal", "Smooth_muscle", "Endothelial", "Endometrial"
# c("Endometrial", "Endothelial", "Fibroblast", "Stromal", "Villious_cytotrophoblast", "Syncytiotrophoblast", "NK", "Macrophage", "Erythroblast")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Trophoblast"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Trophoblast.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Stromal"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Stromal.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Fibroblast"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Fibroblast.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Endothelial"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Endothelial.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Smooth_muscle"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Smooth_muscle.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Endometrial"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Endometrial.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Villious_cytotrophoblast"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Villious_cytotrophoblast.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Syncytiotrophoblast"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Syncytiotrophoblast.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Macrophage"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Macrophage.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("NK"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "NK.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Type"
clusters.0.uninfected <- subset(seurat.object, idents = c("Erythroblast"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Erythroblast.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")
setwd("..")


# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')

# devtools::install_github("hadley/tidyverse")

library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("Type.cluster/"))
manifest<-manifest%>%mutate(name=gsub(".cluster.cellcounts.txt","",value))
manifest
		## These commands below are required, change the prefix to match your samples
df<-read.table("Type.cluster/Trophoblast.cluster.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("cluster","Trophoblast")
 df2<-df
df<-df%>%pivot_longer(cols = `Trophoblast`,names_to = "Type",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("Type.cluster/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("cluster",X)
  df<-df%>%pivot_longer(cols = X,names_to = "Type",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of clusters ## there are 12 Types
# df2%>%mutate(perc=percs_by_group(count,group = Type))
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


final<-final%>%
  mutate(Type=gsub(".cluster.cellcounts.txt","",Type))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = Type))
#  pivot_wider(names_from = Type,values_from = count)

write.table(final,"Type.CellCounts_merged.tsv",sep = "\t",row.names = F)

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

final$cluster<-as.factor(final$cluster)

final$Type<-as.factor(final$Type)
​
final
​
library(ggpubr)
ggbarplot(data = final,
          x = "Type",
          y = "Percent",
          color = "cluster",
          fill = "cluster",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)

png(file=paste0("Type.cluster_diversity.plot.png"),
                res=300, 
                width=1500, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "Type",
          y = "Percent",
          color = "cluster",
          fill = "cluster",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()

####################################################################################################

# reference.list <- all_merged.list[c("Mock-enox-1", "Mock-enox-2", "ZIKV-sham", "ZIKV-enox", "Uninfected")]


dir.create("orig.ident.Type")
setwd("orig.ident.Type")
Idents(seurat.object)<-"orig.ident"
seurat.object2 <- subset(seurat.object, idents = c("Mock-enox-1"))
Idents(seurat.object2)<-"Type"
unfiltered.count <- table(Idents(seurat.object2))
write.table(unfiltered.count, "Mock-enox-1.Type.cellcounts.txt", sep="\t")
Idents(seurat.object)<-"orig.ident"
seurat.object2 <- subset(seurat.object, idents = c("Mock-enox-2"))
Idents(seurat.object2)<-"Type"
unfiltered.count <- table(Idents(seurat.object2))
write.table(unfiltered.count, "Mock-enox-2.Type.cellcounts.txt", sep="\t")
Idents(seurat.object)<-"orig.ident"
seurat.object2 <- subset(seurat.object, idents = c("ZIKV-sham"))
Idents(seurat.object2)<-"Type"
unfiltered.count <- table(Idents(seurat.object2))
write.table(unfiltered.count, "ZIKV-sham.Type.cellcounts.txt", sep="\t")
Idents(seurat.object)<-"orig.ident"
seurat.object2 <- subset(seurat.object, idents = c("ZIKV-enox"))
Idents(seurat.object2)<-"Type"
unfiltered.count <- table(Idents(seurat.object2))
write.table(unfiltered.count, "ZIKV-enox.Type.cellcounts.txt", sep="\t")
Idents(seurat.object)<-"orig.ident"
seurat.object2 <- subset(seurat.object, idents = c("Uninfected"))
Idents(seurat.object2)<-"Type"
unfiltered.count <- table(Idents(seurat.object2))
write.table(unfiltered.count, "Uninfected.Type.cellcounts.txt", sep="\t")
rm(seurat.object2)
setwd("..")


library(tidyverse)
library(mosaic)
manifest<-as_tibble(list.files("orig.ident.Type/"))
manifest<-manifest%>%mutate(name=gsub(".cluster.cellcounts.txt","",value))
manifest

## manually change file names to make sure it works before running function

df<-read.table("orig.ident.Type/Mock-enox-1.Type.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type","AF-p1")
df2<-df
df<-df%>%pivot_longer(cols = `Mock-enox-1`,names_to = "orig.ident",values_to = "Count")
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

## this is required too
df2<-df

## make sure this matches the number of clusters ## there are 22 SubTypes
# df2%>%mutate(perc=percs_by_group(count,group = SubType))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])


​
final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)



​final
final<-final%>%mutate(orig.ident=gsub(".Type.cellcounts.txt","",orig.ident))
final<-final%>%mutate(Percent=percs_by_group(Count,group = orig.ident))
#  pivot_wider(names_from = SubType,values_from = count)
​
write.table(final,"orig.ident.Type.CellCounts_merged.tsv",sep = "\t",row.names = F)

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

## make sure you have enough colors: here, 22
mypal<-get_palette("ucscgb",23)


# final$cluster<-as.factor(final$cluster)
final$orig.ident<-as.factor(final$orig.ident)
final$Type<-as.factor(final$Type)
​
final
​
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
                width=3000, 
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



dir.create("orig.ident.cluster")
setwd("orig.ident.cluster")
# reference.list <- all_merged.list[c("Mock-enox-1", "Mock-enox-2", "ZIKV-sham", "ZIKV-enox", "Uninfected")]


Idents(seurat.object)<-"orig.ident"
clusters.0.uninfected <- subset(seurat.object, idents = c("Mock-enox-1"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.unfilcount <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Mock-enox-1.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"orig.ident"
clusters.0.uninfected <- subset(seurat.object, idents = c("Mock-enox-2"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Mock-enox-2.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"orig.ident"
clusters.0.uninfected <- subset(seurat.object, idents = c("ZIKV-sham"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "ZIKV-sham.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"orig.ident"
clusters.0.uninfected <- subset(seurat.object, idents = c("ZIKV-enox"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "ZIKV-enox.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"orig.ident"
clusters.0.uninfected <- subset(seurat.object, idents = c("Uninfected"))
Idents(clusters.0.uninfected) <- "seurat_clusters"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Uninfected.cluster.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

setwd("..")


library(tidyverse)
library(tidyr)

library(mosaic)
manifest<-as_tibble(list.files("orig.ident.cluster/"))
manifest<-manifest%>%mutate(name=gsub(".cluster.cellcounts.txt","",value))
manifest
        ## These commands below are required, change the prefix to match your samples
df<-read.table("orig.ident.cluster/Mock-enox-1.cluster.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("cluster","Mock-enox-1")
 df2<-df
df<-df%>%pivot_longer(cols = `Mock-enox-1`,names_to = "orig.ident",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("orig.ident.cluster/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("cluster",X)
  df<-df%>%pivot_longer(cols = X,names_to = "orig.ident",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of clusters ## there are 12 orig.idents
# df2%>%mutate(perc=percs_by_group(count,group = orig.ident))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])
final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)

final<-final%>%
  mutate(orig.ident=gsub(".cluster.cellcounts.txt","",orig.ident))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = orig.ident))
#  pivot_wider(names_from = orig.ident,values_from = count)

write.table(final,"orig.ident.CellCounts_merged.tsv",sep = "\t",row.names = F)

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

final$cluster<-as.factor(final$cluster)

final$orig.ident<-as.factor(final$orig.ident)
​
final
​
library(ggpubr)
ggbarplot(data = final,
          x = "orig.ident",
          y = "Percent",
          color = "cluster",
          fill = "cluster",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)

png(file=paste0("orig.ident.cluster_diversity.plot.png"),
                res=300, 
                width=1500, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "orig.ident",
          y = "Percent",
          color = "cluster",
          fill = "cluster",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()





dir.create("seurat_clusters.Type")
setwd("seurat_clusters.Type")

# "Trophoblast", "Fibroblast", "Stromal", "Smooth_muscle", "Endothelial", "Endometrial"
# c("Endometrial", "Endothelial", "Fibroblast", "Stromal", "Villious_cytotrophoblast", "Syncytiotrophoblast", "NK", "Macrophage", "Erythroblast")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("0"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "0.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("1"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "1.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("2"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "2.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("3"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "3.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("4"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "4.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("5"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "5.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("6"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "6.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("7"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "7.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("8"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "8.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("9"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "9.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("10"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "10.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("11"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "11.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("12"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "12.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("13"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "13.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("14"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "14.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("15"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "15.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")


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
df<-read.table("seurat_clusters.Type/1.Type.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type","1")
 df2<-df
df<-df%>%pivot_longer(cols = `1`,names_to = "seurat_clusters",values_to = "Count")
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
## make sure this matches the number of Types ## there are 12 seurat_clusterss
# df2%>%mutate(perc=percs_by_group(count,group = seurat_clusters))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])
df7<-my_fxn(X = manifest$value[6])


??

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
  mutate(seurat_clusters=gsub(".Type.cellcounts.txt","",seurat_clusters))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = seurat_clusters))
#  pivot_wider(names_from = seurat_clusters,values_from = count)

write.table(final,"seurat_clusters.CellCounts_merged.tsv",sep = "\t",row.names = F)

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
​
final
​
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
                width=2500, 
                height=1500)
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

####################################################################################################

####### The clusters are not partitioning clearly by type. Try to increase the resolution and see if that helps?? 

dir.create("res.up")
setwd("res.up")

# seurat.object <- ScaleData(seurat.object, features = integrated.genes, assay = "integrated", vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"))
seurat.object<- RunPCA(seurat.object, features = integrated.genes, assay = "integrated", slot = "scale.data")

seurat.object <- FindNeighbors(seurat.object, features = "integrated.genes", dims = 1:30)
seurat.object <- FindClusters(seurat.object, resolution = 1.7)
seurat.object <- RunUMAP(seurat.object, dims=1:30)
	## with res 1.5 I got 25 clusters and cluster 21 was macrophages with some potentially mis-assigned spatial spots.
p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident")
p1
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE")
p2
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type")
p3
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Subtype")
p4
umap.combined <- CombinePlots(plots = list(p1, p2, p3, p4))
#ggsave("UMAP.pdf", plot = umap.combined, device = "pdf", width = 15, height = 12, units = "in")
ggsave("integrated_UMAP-combined.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("integrated_UMAP-combined-wide.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")

dir.create("seurat_clusters.Type")
setwd("seurat_clusters.Type")

# "Trophoblast", "Fibroblast", "Stromal", "Smooth_muscle", "Endothelial", "Endometrial"
# c("Endometrial", "Endothelial", "Fibroblast", "Stromal", "Villious_cytotrophoblast", "Syncytiotrophoblast", "NK", "Macrophage", "Erythroblast")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("0"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "0.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("1"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "1.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("2"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "2.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("3"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "3.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("4"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "4.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("5"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "5.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("6"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "6.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("7"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "7.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("8"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "8.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("9"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "9.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("10"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "10.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("11"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "11.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("12"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "12.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("13"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "13.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("14"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "14.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("15"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "15.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("16"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "16.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("17"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "17.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("18"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "18.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("19"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "19.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")


Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("20"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "20.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")


Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("21"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "21.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("22"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "22.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("23"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "23.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("24"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "24.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("25"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "25.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")

Idents(seurat.object)<-"seurat_clusters"
Types.0.uninfected <- subset(seurat.object, idents = c("26"))
Idents(Types.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(Types.0.uninfected))
write.table(unfiltered.count, "26.Type.cellcounts.txt", sep="\t")
rm("Types.0.uninfected")



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
df<-read.table("seurat_clusters.Type/1.Type.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type","1")
 df2<-df
df<-df%>%pivot_longer(cols = `1`,names_to = "seurat_clusters",values_to = "Count")
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
## make sure this matches the number of Types ## there are 12 seurat_clusterss
# df2%>%mutate(perc=percs_by_group(count,group = seurat_clusters))
df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])
df7<-my_fxn(X = manifest$value[6])


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


df11<-my_fxn(X = manifest$value[10])
df12<-my_fxn(X = manifest$value[11])
df13<-my_fxn(X = manifest$value[12])
df14<-my_fxn(X = manifest$value[13])
df15<-my_fxn(X = manifest$value[14])
df16<-my_fxn(X = manifest$value[15])
df17<-my_fxn(X = manifest$value[16])
df18<-my_fxn(X = manifest$value[17])
df19<-my_fxn(X = manifest$value[18])
df20<-my_fxn(X = manifest$value[19])
df21<-my_fxn(X = manifest$value[20])
df22<-my_fxn(X = manifest$value[21])
df23<-my_fxn(X = manifest$value[22])
df24<-my_fxn(X = manifest$value[23])
df25<-my_fxn(X = manifest$value[24])
df26<-my_fxn(X = manifest$value[25])
df27<-my_fxn(X = manifest$value[26])
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
final<-full_join(final,df18)
final<-full_join(final,df19)
final<-full_join(final,df20)
final<-full_join(final,df21)
final<-full_join(final,df22)
final<-full_join(final,df23)
final<-full_join(final,df24)
final<-full_join(final,df25)
final<-full_join(final,df26)
final<-full_join(final,df27)

final<-final%>%
  mutate(seurat_clusters=gsub(".Type.cellcounts.txt","",seurat_clusters))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = seurat_clusters))
#  pivot_wider(names_from = seurat_clusters,values_from = count)

write.table(final,"seurat_clusters.CellCounts_merged.tsv",sep = "\t",row.names = F)

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
​
final
​
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
                width=2500, 
                height=1500)
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

####################################################################################################



dir.create("cell.counts")
setwd("cell.counts")

Idents(seurat.object) <- "Condition"
levels(seurat.object)


dir.create("Condition.Type")
setwd("Condition.Type")

Idents(seurat.object)<-"Condition"
clusters.0.uninfected <- subset(seurat.object, idents = c("Control"))
Idents(clusters.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "Uninfected.Type.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Condition"
clusters.0.uninfected <- subset(seurat.object, idents = c("Preeclampsia"))
Idents(clusters.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "ZIKV.Type.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Condition"
clusters.0.uninfected <- subset(seurat.object, idents = c("GDM"))
Idents(clusters.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "GDM.Type.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Condition"
clusters.0.uninfected <- subset(seurat.object, idents = c("SARS-CoV-2"))
Idents(clusters.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "SARS-CoV-2.Type.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Condition"
clusters.0.uninfected <- subset(seurat.object, idents = c("SARS-CoV-2.ax"))
Idents(clusters.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "SARS-CoV-2.Type.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

Idents(seurat.object)<-"Condition"
clusters.0.uninfected <- subset(seurat.object, idents = c("SARS-CoV-2.sx"))
Idents(clusters.0.uninfected) <- "Type"
unfiltered.count <- table(Idents(clusters.0.uninfected))
write.table(unfiltered.count, "SARS-CoV-2.Type.cellcounts.txt", sep="\t")
rm("clusters.0.uninfected")

setwd("..")

# install.packages('tidyverse')
# install.packages('mosaic')
# install.packages('ggpubr')
# devtools::install_github("hadley/tidyverse")
library(tidyverse)
library(tidyr)
library(mosaic)
manifest<-as_tibble(list.files("Condition.Type/"))
manifest<-manifest%>%mutate(name=gsub(".Type.cellcounts.txt","",value))
manifest
      ## These commands below are required for at least 1 of the data frames. Make sure to change the file pre-fix if necessary [AF_Epithelial]
df<-read.table("Condition.Type/Control.Type.cellcounts.txt",header = T,sep = "\t",row.names = 1)
colnames(df)<-c("Type","Control")
 df2<-df
df<-df%>%pivot_longer(cols = `Control`,names_to = "Condition",values_to = "Count")
df
manifest

my_fxn<-function(X){
  df<-read.table(paste0("Condition.Type/",X),
                 header = T,
                 sep = "\t",
                 row.names = 1)
  colnames(df)<-c("Type",X)
  df<-df%>%pivot_longer(cols = X,names_to = "Condition",values_to = "Count")
  df2<-full_join(df,df2)
  return(df)
}

df2<-df
## make sure this matches the number of groups

df2<-my_fxn(X = manifest$value[1])
df3<-my_fxn(X = manifest$value[2])
df4<-my_fxn(X = manifest$value[3])
df5<-my_fxn(X = manifest$value[4])
df6<-my_fxn(X = manifest$value[5])
df7<-my_fxn(X = manifest$value[6])

final<-full_join(df2,df3)
final<-full_join(final,df4)
final<-full_join(final,df5)
final<-full_join(final,df6)
final<-full_join(final,df7)


final<-final%>%
  mutate(SampleID=gsub(".Type.cellcounts.txt","",SampleID))
final<-final%>%
  mutate(Percent=percs_by_group(Count,group = SampleID))
write.table(final,"Condition.Type.CellCounts_merged.tsv",sep = "\t",row.names = F)
final
library(ggpubr)
mypal<-get_palette("igv",22)
final$Type<-as.factor(final$Type)
final$SampleID<-as.factor(final$SampleID)
final
library(ggpubr)
ggbarplot(data = final,
          x = "Condition",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)

png(file=paste0("Condition.Type_diversity.plot.png"),
                res=300, 
                width=2000, 
                height=2500)
feature.plot <- ggbarplot(data = final,
          x = "Condition",
          y = "Percent",
          color = "Type",
          fill = "Type",
          position = position_stack(),
          ggtheme = theme_pubr(),
          palette=mypal)
feature.plot + theme(axis.text.x = element_text(angle = 90))
dev.off()

