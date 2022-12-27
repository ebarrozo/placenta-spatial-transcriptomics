## SpaceRanger_ERB_v1.sh
# Enrico Barrozo, BCM, Aagaard Lab, 2021

## Analyze data on AagaardLab3
# 


## download latest spaceranger and mouse reference
## https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/installation
curl -o spaceranger-1.3.0.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-1.3.0.tar.gz?Expires=1634802546&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvc3BhdGlhbC1leHAvc3BhY2VyYW5nZXItMS4zLjAudGFyLmd6IiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNjM0ODAyNTQ2fX19XX0_&Signature=TFCSeNH0mg4Ym9ugg3dC~a41SRlGCP--XLaxFc~1tXFYA-TMaGVZgsMT1OKS9RLlxrkftF0IVSVPID4VElOt0hdapliDLGIGZazBueFfyWAlElOWEmkJdIvJUBKlN8cB04YQ9Gg~73Dt8LdPMUdPb4IFGKqG6rlZvOeJDqwKrKuZ1Yfo0GXC~eoIVXwNdSU9CBI7DGhyE1VGEOff6Hl9YVh~crYLFRTrZJ808aLD0UuT~HkdULp1P9dsrnxdmXUBNaKdmCBw5sBqyiB95adw5Gnyxyh-a2zmY6qlp01vlOcI0R4xs6lb4dTnq-Inxpr2sgldZHt1Q6pD3QiwBR-7uw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

curl -O https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-mm10-2020-A.tar.gz

cd
tar -xzvf spaceranger-1.3.0.tar.gz
tar -xzvf refdata-gex-mm10-2020-A.tar.gz

	## human + SARS-CoV-2 can now  be found at /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2
################################################# Alignment in SpaceRanger ###########################################################
cd
export PATH=/home/ebarrozo/cellranger-6.1.1:$PATH
export PATH=/home/ebarrozo/spaceranger-1.3.0:$PATH

cd /home/ebarrozo/visium/docs
	## make sure the .tif files of each slide from the pathology core are transfered here

cd /home/ebarrozo/visium/data
	## make sure are all fastq files are named correctly
			# Consistent with bcl2fastq/mkfastq, e.g. "MySample_S1_L001_R1_001.fastq.gz" "MySample_S1_L001_R2_001.fastq.gz"


cd /home/ebarrozo/visium/results
## mkfastq done by GARP/multiomics core
	# spaceranger mkfastq

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
cd /home/ebarrozo/visium/results
spaceranger count --id=S01-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25749_01_HCV \
                   --sample=KA_25749 \
                   --image=/home/ebarrozo/visium/docs/1.tif \
                   --slide=V10S21-369 \
                   --area=A1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-369-A1.json \
                   --localcores=48 \
                   --localmem=300

## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)

################################################# ################################################# 
################################################# I found after running auto and manual with S07 that manual is essential. Auto image alignment was traaash. 
################################################# 
spaceranger count --id=S01-auto \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25749_01_HCV \
                   --sample=KA_25749 \
                   --image=/home/ebarrozo/visium/docs/1.tif \
                   --slide=V10S21-369 \
                   --area=A1 \
                   --localcores=48 \
                   --localmem=300
# Downloading a Slide File for Local Operation; https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count
################################################# 

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S03-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25750_03_HD \
                   --sample=KA_25750 \
                   --image=/home/ebarrozo/visium/docs/3.tif \
                   --slide=V10S21-369 \
                   --area=B1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-369-B1.json \
                   --localcores=48 \
                   --localmem=300

## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S03-auto \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25750_03_HD \
                   --sample=KA_25750 \
                   --image=/home/ebarrozo/visium/docs/3.tif \
                   --slide=V10S21-369 \
                   --area=B1 \
                   --localcores=48 \
                   --localmem=300
################################################# 

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S04-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25751_04_HCMR \
                   --sample=KA_25751 \
                   --image=/home/ebarrozo/visium/docs/4.tif \
                   --slide=V10S21-369 \
                   --area=C1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-369-C1.json \
                   --localcores=48 \
                   --localmem=300

## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S04-auto \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25751_04_HCMR \
                   --sample=KA_25751 \
                   --image=/home/ebarrozo/visium/docs/4.tif \
                   --slide=V10S21-369 \
                   --area=C1 \
                   --localcores=48 \
                   --localmem=300                   
################################################# 

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S15-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25756_15_HCV \
                   --sample=KA_25756 \
                   --image=/home/ebarrozo/visium/docs/15.tif \
                   --slide=V10S21-368 \
                   --area=D1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-368-D1.json \
                   --localcores=48 \
                   --localmem=300

## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S15-auto \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25756_15_HCV \
                   --sample=KA_25756 \
                   --image=/home/ebarrozo/visium/docs/15.tif \
                   --slide=V10S21-368 \
                   --area=D1 \
                   --localcores=48 \
                   --localmem=300
################################################# 



      ## make sure are all fastq files are named correctly
                  # Consistent with bcl2fastq/mkfastq, e.g. "MySample_S1_L001_R1_001.fastq.gz" "MySample_S1_L001_R2_001.fastq.gz"

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S16-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo//aagaard-tillery_134233 \
                   --sample=KA26045 \
                   --image=/home/ebarrozo/visium/docs/16.tif \
                   --slide=V10S21-370 \
                   --area=A1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-370-A1.json \
                   --localcores=48 \
                   --localmem=300

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S17-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo//aagaard-tillery_134233 \
                   --sample=KA26046 \
                   --image=/home/ebarrozo/visium/docs/17.tif \
                   --slide=V10S21-370 \
                   --area=B1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-370-B1.json \
                   --localcores=48 \
                   --localmem=300

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S18-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo//aagaard-tillery_134233 \
                   --sample=KA26047 \
                   --image=/home/ebarrozo/visium/docs/18.tif \
                   --slide=V10S21-370 \
                   --area=C1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-370-C1.json \
                   --localcores=48 \
                   --localmem=300

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S19-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo//aagaard-tillery_134233 \
                   --sample=KA26048 \
                   --image=/home/ebarrozo/visium/docs/19.tif \
                   --slide=V10S21-370 \
                   --area=D1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-370-D1.json \
                   --localcores=48 \
                   --localmem=300

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S21-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo//aagaard-tillery_134233 \
                   --sample=KA26049 \
                   --image=/home/ebarrozo/visium/docs/21.tif \
                   --slide=V10S21-371 \
                   --area=A1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-371-A1.json \
                   --localcores=48 \
                   --localmem=300

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S22-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo//aagaard-tillery_134233 \
                   --sample=KA26050 \
                   --image=/home/ebarrozo/visium/docs/22.tif \
                   --slide=V10S21-371 \
                   --area=B1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-371-B1.json \
                   --localcores=48 \
                   --localmem=300

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S23a-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo//aagaard-tillery_134233 \
                   --sample=KA26051 \
                   --image=/home/ebarrozo/visium/docs/23a.tif \
                   --slide=V10S21-371 \
                   --area=C1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-371-C1.json \
                   --localcores=48 \
                   --localmem=300

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S23b-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo//aagaard-tillery_134233 \
                   --sample=KA26052 \
                   --image=/home/ebarrozo/visium/docs/23b.tif \
                   --slide=V10S21-371 \
                   --area=D1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-371-D1.json \
                   --localcores=48 \
                   --localmem=300








###################################################### GE Slide 5, change fastq folder and sample names. 


################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S20-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo//aagaard-tillery_134233 \
                   --sample=KA_25756 \
                   --image=/home/ebarrozo/visium/docs/20.tif \
                   --slide=V10N16-022 \
                   --area=A1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10N16-022-A1.json \
                   --localcores=48 \
                   --localmem=300
################################################# 

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S24-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_134233 \
                   --sample=KA_25756 \
                   --image=/home/ebarrozo/visium/docs/24.tif \
                   --slide=V10N16-022 \
                   --area=B1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10S21-368-B1.json \
                   --localcores=48 \
                   --localmem=300
################################################# 

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S25-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25756_25_HCV \
                   --sample=KA_25756 \
                   --image=/home/ebarrozo/visium/docs/25.tif \
                   --slide=V10N16-022 \
                   --area=C1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10N16-022-C1.json \
                   --localcores=48 \
                   --localmem=300
################################################# 

################################################# 
## Run w/ and w/o manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S26-manual \
                   --transcriptome=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2 \
                   --fastqs=/media/jochum00/Aagaard_Raid2/ebarrozo/aagaard-tillery_133875/25756_26_HCV \
                   --sample=KA_25756 \
                   --image=/home/ebarrozo/visium/docs/26.tif \
                   --slide=V10N16-022 \
                   --area=D1 \
                   --loupe-alignment=/home/ebarrozo/visium/docs/manual_image_alignments/V10N16-022-D1.json \
                   --localcores=48 \
                   --localmem=300
################################################# 

A successful spaceranger count run concludes with a message similar to this:
2016-11-10 16:10:09 [runtime] (join_complete)   ID.sample345.SPATIAL_RNA_COUNTER_CS.SPATIAL_RNA_COUNTER_CS.SUMMARIZE_REPORTS
 
Outputs:
- Run summary HTML:                         /opt/sample345/outs/web_summary.html
- Outputs of spatial pipeline:              /opt/sample345/outs/spatial
- Run summary CSV:                          /opt/sample345/outs/metrics_summary.csv
- BAM:                                      /opt/sample345/outs/possorted_genome_bam.bam
- BAM index:                                /opt/sample345/outs/possorted_genome_bam.bam.bai
- Filtered feature-barcode matrices MEX:    /opt/sample345/outs/filtered_feature_bc_matrix
- Filtered feature-barcode matrices HDF5:   /opt/sample345/outs/filtered_feature_bc_matrix.h5
- Unfiltered feature-barcode matrices MEX:  /opt/sample345/outs/raw_feature_bc_matrix
- Unfiltered feature-barcode matrices HDF5: /opt/sample345/outs/raw_feature_bc_matrix.h5
- Secondary analysis output CSV:            /opt/sample345/outs/analysis
- Per-molecule read information:            /opt/sample345/outs/molecule_info.h5
- Loupe Browser file:                       /opt/sample345/outs/cloupe.cloupe
# - Spatial Enrichment using Moran's I file:  /opt/sample345/outs/spatial_enrichment.csv
 
Pipestance completed successfully!
