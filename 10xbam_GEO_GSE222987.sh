# 10xbam_GEO.sh

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-FPlatform Medicine

# GEO Contents: Spatial transcriptomics data are available at GEO Accession GSE222987. Use the 10x Genomics bam2fastq package for manual re-alignments. Processed .h5 matricies are also deposited. The script used to generate the custom human + SARS-CoV-2 reference (GRCh38_and_NC_045512.2) is provided here. For each sample, we have uploaded the 10x BAM file (possorted_genome_bam.bam) and processed count matrix (filtered_feature_bc_matrix.h5). Do not use fastq-dump or sam-dump to pull the data as the 10x Bam is different frum regular BAM files. The 10x BAM is different from regular BAM files, so the sam-dump option does not work. Please download the BAM file with the link in the 'Data access' --> 'Original format' section (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR23094132&display=data-access). Please see the 10x Genomics Support Page on the topic: https://kb.10xgenomics.com/hc/en-us/articles/360035250172-Why-do-I-not-see-any-FASTQs-or-incomplete-FASTQs-for-the-SRA-of-my-interest- . Then run bamtofastq (https://support.10xgenomics.com/docs/bamtofastq) in the CellRanger package to extract all the fastq files in the format necessary for CellRanger alignments. For example, 1 10x BAM produces 12 total fastq files. 
## for example- first link is where to find the 'original file' download link, and on the second is the link for download.
   # https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR23094132&display=data-access	
   # https://sra-pub-src-1.s3.amazonaws.com/SRR23094133/S25_possorted_genome_bam.bam.1

#For example- the 10x BAM for 1 sample generates 12 fastq files, in the format necessary for CellRanger alignment. 
  # Therefore, we uploaded the 16 10x BAM files instead of the 192 fastq files. 

## These are all of the SRR samples 
# SRR23094132 SRR23094133 SRR23094134 SRR23094135 SRR23094136 SRR23094137 SRR23094138 SRR23094139 SRR23094140 SRR23094141 SRR23094142 SRR23094143 SRR23094144 SRR23094145 SRR23094146 SRR23094147 

## pull all of the original files
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094133/S25_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094133/S25_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094134/S24_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094135/S23b_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094136/S23a_possorted_genome_bam.bam.1 
wget ttps://sra-pub-src-1.s3.amazonaws.com/SRR23094137/S22_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094138/S21_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094139/S20_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094140/S19_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094141/S18_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094142/S17_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094143/S16_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094144/S15_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094145/S04_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094146/S03_possorted_genome_bam.bam.1 
wget https://sra-pub-src-1.s3.amazonaws.com/SRR23094147/S01_possorted_genome_bam.bam.1 

## convert the bam files to fastqs
for f in $.bam.1; do
	bamtofastq $f . 
done




