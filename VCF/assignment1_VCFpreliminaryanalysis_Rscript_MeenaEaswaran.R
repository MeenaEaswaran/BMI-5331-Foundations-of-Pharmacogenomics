#set working directory
#setwd()

#note the system time 
Sys.time()

#install vcfR library if required
#install.packages("vcfR")

#load vcfR library 
library(vcfR)

#read the genotype VCF file 
vcf <- read.vcfR("assignment1_geno.vcf", verbose = FALSE)
vcf

# Check the first 6 lines of the VCF to get information on the meta, fixed and genotype sections
head(vcf)

#Preliminary analysis with VCF:
#a.	Does the file contain diploid alleles or haploid?
#Based on the above commands, the file seems to contain diploid alleles (2 chromosome sets)
#In a VCF file, diploid genotypes are typically represented as two alleles separated by either “/” or “|”
#Confirm with the commands below

# Create a data frame that maps the unique IDs to the original IDs 
#unique ID needed for genotype extraction as ID contains non-unique names (duplicates)
id_map <- data.frame(
  UniqueID = paste0("var", seq_len(nrow(vcf@fix))),
  OriginalID = vcf@fix[,"ID"]
)
#Unique IDS Will be needed for mapping original variants at step 3

#Check for potential duplicates in the original variant IDs
duplicates_origignalvariantID <- id_map[duplicated(id_map$OriginalID) | duplicated(id_map$OriginalID, fromLast = TRUE),]
duplicates_origignalvariantID
#duplicates present for one variant rs112607901

# Assign unique IDs to the variants 
vcf@fix[,"ID"] <- id_map$UniqueID

# Extract genotype data with element GT
geno_data <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

# Check if the genotypes are diploid or haploid
if (any(grepl("/|\\|", geno_data))) {
  print("The VCF file contains diploid alleles.")
} else {
  print("The VCF file contains haploid alleles.")
}

#Alternatively, check for unique genotypes
unique(as.vector(geno_data)) 
#The outcomes give unique genotypes in terms of ploidy in the VCF file 
#all are diploid and phased

# Convert geno_data to a data frame 
geno_data <- as.data.frame(geno_data)

#################################################################

#b. To check if the alleles are phased or unphased:

#Phased alleles are denoted by | and Unphased alleles are denoted by /
#Based on the above command, when checked for unique genotypes, all alleles are phased
#Confirm by using the following commands

# Check if the genotype data contains the phasing character "|"
phased_info <- grep("|", vcf@gt, value = TRUE)

# If the length of the phased_info vector is greater than 0, the alleles are phased
if(length(phased_info) > 0) {
  cat("The alleles are phased.")
} else {
  cat("The alleles are unphased.")
}

#################################################################

#c. To count how many variants have passed quality control:

# Check unique quality scores and unique filter notations 
unique_qual_scores <- unique(vcf@fix[, "QUAL"])
unique_qual_scores

unique_filter <- unique(vcf@fix[, "FILTER"])
unique_filter

#Based on above command, 100 is the only unique QUAL score and PASS is the only unique FILTER notation
#This preliminary indicates all 6500 variants passed the quality control. Rule out the role of NAs as below

# Count the number of variants that have a non-NA quality score
passed_variants <- sum(!is.na(vcf@fix[, "QUAL"]))
passed_variants

#################################################################

#d. To count how many SNPs exist in the file

vcf2 <- extract.indels(vcf, return.indels = FALSE) 
#Setting return.indels false, will return only SNPs. If TRUE, it will return only indels
vcf2

#Extract information
info_vcf2 <- extract_info_tidy(vcf2)

#Ensure if SNP is only unique form in the variant type column
unique(info_vcf2$VT)

#to know the length or count of the revised VCF
SNP_count <- length(vcf2@fix[, "ID"])
SNP_count

#################################################################

#e.	How many indels exist in the file? Can you tell how many of them are deletions?

vcf3 <- extract.indels(vcf, return.indels = TRUE) 
#Setting return.indels TRUE, will return only indels
vcf3

#Extract information 
info_vcf3<- extract_info_tidy(vcf3)

#Ensure if indel is only unique form in the variant type column
#This should have no SNP
unique(info_vcf3$VT)

#this is not the case as "INDEL", "SNP,INDEL" "SV"  do come up

#to get only indel count
# Get unique counts for each type of variant present (, "INDEL", "SNP,INDEL", "SV")
counts_vcf3 <- table(info_vcf3$VT)
counts_vcf3

#explore metadata
vcf@meta
#line 249,  have annotations in Ancestral Allele (AA) column about indels; specific insertions or deletions present

# Install dplyr if required
#install.packages("dplyr")

# Load dplyr package
library(dplyr)

# Extract items with deletion mentioned in the AA column
deletion_rows <- info_vcf3 %>%
  filter(grepl("deletion", AA, ignore.case = TRUE))
nrow(deletion_rows)

#confirm if deletion_rows are part of only indel category and not SNP, INDEL
unique(deletion_rows$VT)

#################################################################

#e.For structural variations (copy number variations)

#explore metadata
vcf@meta
#VCF file has copy number variations as seen in metadata

#CNV information to be obtained from INFO section
info_vcf <- extract_info_tidy(vcf)
unique(info_vcf$VT) 
#confirms presence of structural variants (SV)

# Get unique counts for each type of variant present ("SNP", "INDEL", "SNP,INDEL", "SV")
counts_vcf <- table(info_vcf$VT)
counts_vcf

# Get the count specifically for "SV"
sv_count <- counts_vcf["SV"]
sv_count

# Alternatively confirm  number of structural variant type DEL/DUP associated to CNV in file
#Deletions and duplications with a size larger than 1,000 bases (>1 kb) are often referred to as copy-number variations (CNVs)
c <- table(info_vcf$SVTYPE)
c

#################################################################

#g. To count variants with more than one alternative allele
# based on the VCF meta data, total number of alternate alleles in called genotypes is in line 239 under AC 

# Install stringr if required
#install.packages("stringr")

# Load stringr packages
library(stringr)

# Filter the 'info' data frame to create 'multi_allele_rows', which includes only rows where the 'AC' column contains more than one alternate allele (indicated by the presence of a comma)
multi_allele_rows <- info_vcf %>%
  filter(str_count(AC, ",") > 0)

nrow(multi_allele_rows)

