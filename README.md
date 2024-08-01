# RNA-seq
# Project Title
A computational approach to identifying druggable targets and repurposing existing compounds for the treatment of CNS leukaemia (MSc Biomedical Sciences Thesis Project)

# Overview
This repository contains scripts and data used for the analysis of RNA sequencing (RNA-Seq) data to identify potential drug repurposing candidates for paediatric acute lymphoblastic leukaemia (ALL) with a focus on central nervous system (CNS) involvement.

# Datasets
In vivo dataset: Human T-ALL cell line CCRF-CEM was engrafted into immunodeficient NSG mice. Leukemic blasts were retrieved from the CNS and BM at clinical endpoint. RNA was extracted, and sequencing was performed on an Illumina NextSeq500 system. This dataset includes raw reads in FASTQ files from four BM and four CNS tumour samples.

Patient Dataset: Public RNA-seq dataset (GSE81518) with gene expression profiles from CSF and BM of three B-ALL patients was analyzed for comparison and validation

# Data analysis and visualisation
All data processing tasks were performed using PuTTY and a High-Performance Computing (HPC) server (MARS) supported by the University of Glasgow. Analysis and visualization were conducted in Python and R.

# Workflow 
## 1. Quality control
The first step of the data processing for mouse data was quality control of the raw reads. Packages like FastQC, MultiQC and Trimmomatic were used to remove low-quality reads and trim adapter sequences. The code is included in the "RNAseq_QC.sh" file.

## 2. Removal of mouse reads
The xengsort package was used in "xengsort.sh" to filter out mouse tissue reads, while retaining human tumour reads.

## 3. Mapping and Alignment
Mapping and alignment to the human genome was conducted with the package HISAT2 and the code is available in the "RNAseq_Mapping.sh" file. 

## 4. Quantification and Expression Matrix
The code for the quantification of reads associated with protein-coding genes and the creation of a gene-level count matrix is included in the files "RNAseq_counts.sh" and "RNAseq_matrix.sh" respectively. 

## 5. Data normalisation and identification of outliers, differential expression analysis, functional enrichment analysis, GSEA and data visualisation
The R script "RNASeq_Rscript.R" contains the code for data normalisation with DESeq2 and edgeR, identification of outliers with the VST package, differential expression analysis, functional enrichment analysis and GSEA, as well as data visualisation. This code was used for both in vivo and patient datasets.

## 6. Drug repurposing with ChEMBL
The "CHEMBL.py" and "CHEMBL_top_compounds.py" files contain the code that was used to identify potential compounds that target the most significant genes of each dataset. 

## 7. Protein-protein interaction (PPI) networks 
PPI networks of the 500 most significant genes were generated by using the STRING online database. The STRING PPI networks were further analysed and visualised in Cytoscape (version 3.10.2) with the CytoHubba plug-in to identify protein hubs based on node degrees.

## 8. Drug repurposing with CMap 
Potential drugs for repurposing were identified using the CMap database from Broad Institute. The top 150 upregulated and 150 downregulated DEGs of each dataset were used to generate a list of compounds ranked by connectivity score. Drugs with a connectivity score ≤ -90 were considered for repurposing. 

# Conclusion
This workflow facilitates the identification of potential gene targets and drug repurposing candidates for treating paediatric CNS ALL. By integrating various bioinformatics tools and datasets, this approach aims to enhance treatment efficacy and patient outcomes. 

