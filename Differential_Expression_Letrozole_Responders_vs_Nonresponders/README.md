# Differential Expression Analysis: Responders vs. Non-Responders to Letrozole Treatment 
This folder contains a R script for the analysis of gene expression between two phenotypes—responders and non-responders to Letrozole drug treatment in the context of breast cancer tumors. 
The study leveraged publicly available microarray dataset to identify **differentially expressed genes (DEGs)** between responders and non-responders.

## Methods and Tools
### Data Retrieval
- Source: NCBI GEO
- Dataset: [GDS3116](https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3116)
- Tool Used: GEOQuery (R package)
       
### Differential Gene Expression Analysis
- Tool Used: limma (Linear Models for Microarray Data)
- Statistical Corrections:
  - Bonferroni correction (α = 0.01)
  - Benjamini-Hochberg FDR correction (α = 0.01)

### Mapping Probes to Offical Gene Names
- Tool Used: **[DAVID](https://david.ncifcrf.gov/conversion.jsp)** Gene ID Conversion Tool 
- Identifier: AFFEMETRIX_3PRIME_IVT_ID
- Species: Homo sapiens

## Citation
If you use the tools or dataset mentioned in this repository in your research, please cite the following references:

- Barrett, T., Wilhite, S. E., Ledoux, P., Evangelista, C., Kim, I. F., Tomashevsky, M., Marshall, K. A., Phillippy, K. H., Sherman, P. M., Holko, M., Yefanov, A., Lee, H., Zhang, N., Robertson, C. L., Serova, N., Davis, S., & Soboleva, A. (2013). NCBI GEO: archive for functional genomics data sets—update. Nucleic Acids Research, 41(D1), D991–D995. https://doi.org/10.1093/NAR/GKS1193

- Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res [Internet]. 2015 Apr 20 [cited 2023 Dec 7];43(7):e47–e47. Available from: https://dx.doi.org/10.1093/nar/gkv007

- Sean D, Meltzer PS. GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. Bioinformatics [Internet]. 2007 Jul 15 [cited 2023 Dec 7];23(14):1846–7. Available from: https://dx.doi.org/10.1093/bioinformatics/btm254

- Sherman, B. T., Hao, M., Qiu, J., Jiao, X., Baseler, M. W., Lane, H. C., Imamichi, T., & Chang, W. (2022). DAVID: a web server for functional enrichment analysis and functional annotation of gene lists (2021 update). Nucleic acids research, 50(W1), W216–W221. https://doi.org/10.1093/nar/gkac194

