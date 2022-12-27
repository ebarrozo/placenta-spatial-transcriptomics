# Human placenta spatial transcriptomics analysis

Spatial transcriptomics. Human placentae from distinct regions including the chorionic villi, decidua, and chorioamniotic membranes, or cross-sections from the parenchyma, were fresh-frozen in optimal cutting temperature solution (FF-OCT). Blocks were cryosectioned and H&E stained directly on 10X Genomics Visium Gene Expression slides (v2) and imaged using a Nikon Eclipse SE Ni microscope using a Nikon DS-Ri1 camera, Nikon Plan Apo objective at 10x magnification (0.45 aperture, 0.91 μm/pixel resolution). Tissues were permeabilized and RNA was subject to spatial transcriptomics library preparation including poly(dT) reverse transcription. Libraries were sequenced on the Illumina NovaSeq S4 platform with 2% PhiX. 
	
Single-cell and spatial transcriptomics analyses. Reads were demultiplexed and aligned to a custom human (GRCh38) and SARS-CoV-2 reference genome (NC_045512.2) using SpaceRanger (v1.3.0) and custom bash scripts. Downstream analyses were done using the package Seurat (v4.0.3) in rStudio (v4.1.1). Counts matrices were filtered iteratively to exclude low-quality transcriptomes and clusters defined by quality control metrics (e.g. mitochondrial or hemoglobin gene expression). Spatial transcriptomes were normalized and scaled using a negative binomial model (SCTransform) and the top 3,000 most variable transcripts were used for to principal component analysis dimension reduction. The first 30 PCA dimensions were used for K-nearest neighbors’ analysis, clustered using a Louvain algorithm with the default resolution parameter (0.6), and visualized by unique manifold approximation and projection (UMAP) in two dimensions. Significantly upregulated transcripts were manually examined and compared to the term placenta atlas in addition to using EnrichR(58), PlacentaCellEnrich(110), and the Human Protein Atlas(111) to annotate clusters. Pseudotime trajectory analysis was done using Monocle3 (v1.0.0) utilizing gene sets and starting points denoted in the text. Spatial and scRNA-seq datasets were integrated using reciprocal PCA of 3,000 reference transcripts prior to clustering.

GEO Contents: Spatial transcriptomics data are available at GEO Accession GSEXXXXXX. Use the 10x Genomics bam2fastq package for manual re-alignments. Processed .h5 matricies are also deposited. The script used to generate the custom human + SARS-CoV-2 reference (GRCh38_and_NC_045512.2) is provided here.

SampleIDs: Case-matched metadata are available with in the Supplementary Information. S01 chorionic villi, S03 decidua, and S04 chorioamniotic membrane roll samples are from the same uninfected placenta. S24, S25, and S26 are uninfected controls. S15, S17, S18, and S19 are from maternal asymptomatic (no COVID-19 symptoms) SARS-CoV-2 positive placentae and S16, S20, S21, S22, and S23 are from symptomatic maternal COVID-19 placentae. S23 was sampled twice (S23a and S23b).

# Repository contents
A) Images necessary for alignment of spatial transcriptomes. Higher-resolution (~100MB per tiff) are availible upon request.

B) Bash scripts for creating a custom reference transcriptome and alignment using Spaceranger

C) R script for QC filtering and initial clustering of spatial transcriptomics data using Seurat and initial analyses.

D) Clustering of previously published placenta scRNA-seq and snRNA-seq data using Seurat, and integration of placenta spatial transcriptomics, scRNA-seq, and snRNA-seq data.

E) R script for further analysis of spatial transcriptomics data

F) Monocle pseudotime trajectory analysis of macrophage niches
