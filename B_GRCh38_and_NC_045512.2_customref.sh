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

################################################# Make Custom References ###########################################################
## Download the latest and greatest references for each virus (.gtf and .fasta)
	## SARS-CoV-2 is NC_045512.2
	## ZIKV clinical strain we are using is NC_045512.2 ; revised gtf annotation details can be found at https://github.com/ebarrozo/Aagaard_ZIKV_AGOHITSCLIP/blob/main/A_zikv-hg38_eb_v3.sh
cd /home/ebarrozo
mkdir visium
cd visium
mkdir data
mkdir docs
mkdir scripts
mkdir results

cd /home/ebarrozo/visium/docs/customref
cp /home/ebarrozo/BALF/docs/customref/NC_045512.2.genes.filtered.gtf .
cp /home/ebarrozo/BALF/docs/customref/NC_045512.2.fasta .

# convert to gtf

gffread NC_045512.2.gff3 -T -o NC_045512.2.gtf
grep "exon" NC_045512.2.gtf > NC_045512.2.genes.gtf
	## only 2 genes with exon included... grep "CDS" NC_001802.1.gtf > NC_001802.1.genes.gtf
	## 13 genes included, overlap in ORF1ab 3 times. 

## replace CDS with exon
cat NC_045512.2.genes.filtered.gtf
sed -i 's/CDS/exon/g' NC_045512.2.genes.filtered.gtf
cat NC_045512.2.genes.filtered.gtf


cd
export PATH=/home/ebarrozo/cellranger-6.1.1:$PATH
cd /home/ebarrozo/visium/docs/customref
cellranger mkgtf NC_045512.2.genes.gtf NC_045512.2.genes.filtered.gtf --attribute=gene_biotype:protein_coding
mkgtf NC_001803.1.genes.gtf NC_001803.1.genes.filtered.gtf --attribute=gene_biotype:protein_coding

## edit NC_001803.1.genes.filtered.gtf so ORF1ab is not repeated and transcripts do not overlap
emacs NC_001803.1.genes.filtered.gtf 
	## to exit emacs: ctrl+x then ctrl+c, press y to save


## Filter gtf files to make sure they are formatted for cellranger. Extracts transcripts with "exon" and "protein_coding".
## download and install cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
cd
export PATH=/home/ebarrozo/cellranger-6.1.1:$PATH
cellranger mkgtf NC_001802.1.genes.gtf NC_001802.1.genes.filtered.gtf --attribute=gene_biotype:protein_coding

cd /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs
#################################################################################### HAVE TO USE SPACERANGER MKREF
## Make custom ref
spaceranger mkref --nthreads=48 --genome=GRCh38 \
--fasta=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
--genes=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
--genome=NC_045512.2 \
--fasta=/home/ebarrozo/visium/docs/customref/NC_045512.2.fasta \
--genes=/home/ebarrozo/visium/docs/customref/NC_045512.2.genes.filtered.gtf
	## human + SARS-CoV-2 can now  be found at /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_045512.2


