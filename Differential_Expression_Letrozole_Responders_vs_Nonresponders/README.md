# Differential Expression Analysis: Responders vs. Non-Responders to Letrozole Treatment 
This folder contains R script for the analysis of gene expression between two phenotypes—responders and non-responders to Letrozole drug treatment in the context of breast cancer tumors. 
The study leverages publicly available data to identify **differentially expressed genes (DEGs)**.

## Methods and Tools
### Data Retrieval
- Source: NCBI GEO
- Dataset: GDS3116 (https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3116)
- Tool Used: GEOQuery (R package)
       
### Differential Gene Expression Analysis
- Tool Used: limma (Linear Models for Microarray Data)
- Statistical Corrections:
  - Bonferroni correction (α = 0.01)
  - Benjamini-Hochberg FDR correction (α = 0.01)

### Mapping Probes to Offical Gene Names
- Tool Used: **DAVID** Gene ID Conversion Tool (https://david.ncifcrf.gov/conversion.jsp)
- Identifier: AFFEMETRIX_3PRIME_IVT_ID
- Species: Homo sapiens
