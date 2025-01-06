# Differential Expression Analysis: Responders vs. Non-Responders to Letrozole Treatment 
This folder contains the analysis of gene expression between two phenotypes—responders and non-responders to Letrozole drug treatment in the context of breast cancer tumors. 
The study leverages publicly available data and bioinformatics tools to identify **differentially expressed genes (DEGs)**.

## Methods and Tools
1. ### Data Retrieval
       Source: NCBI GEO
       Dataset: GDS3116 (https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3116)
       Tool Used: GEOQuery (R package)
       
2. ### Differential Gene Expression Analysis
       Tool Used: limma (Linear Models for Microarray Data)
       Statistical Corrections:
       Bonferroni correction (α = 0.01)
       Benjamini-Hochberg FDR correction (α = 0.01)
