#set working directory 
#setwd()

#note the system time 
Sys.time()

#install GEOquery package if required
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("GEOquery", force = TRUE)

#load GEOquery library
library(GEOquery)

# Download and read the DataSet SOFT file

gds <- getGEO(GEO = "GDS3116", filename = "GDS3116.soft.gz")
#gds <- getGEO("GDS3116", GSEMatrix = FALSE, AnnotGPL = FALSE) #This should work too

############### Preliminary exploration of file to familiarize with GEO dataset #############################

#Explore gds meta data 
Meta(gds)$description
Meta(gds)$platform
Meta(gds)$sample_count
Meta(gds)$sample_organism
Meta(gds)$sample_type
Meta(gds)$title

#Explore sample information associated to each phenotype 
gds@dataTable@columns[["description"]]

#Explore probe information if required
#gds@dataTable@table[["ID_REF"]]

############### Convert GDS to ExpressionSet #############################

#Convert GDS to an ExpressionSet
eset <- GDS2eSet(gds,do.log2=FALSE) #do.log2 as TRUE will log transform

#View phenotype data in Expressionset as a data frame 
phenotype_df <- as.data.frame (pData(eset))
View(phenotype_df)
#This file shows that required information about two relevant phenotypes: "responder" and "nonresponder" are in the description column

############### Perform phenotype filtering based on description column #############################

#1: Convert ExpressionSet as a data frame
eset_df <- as.data.frame(eset) 
View(eset_df)
#This will have X prefixed to all probes IDs in the data frame columns except for the last 4 columns

# Retrieve actual IDs and assign them as column names of the data frame
colnames(eset_df) <- gds@dataTable@table[["ID_REF"]]

# Assign specific names as earlier in the phenotype data to the last four columns
colnames(eset_df)[(ncol(eset_df) - 3):ncol(eset_df)] <- c("sample", "agent", "individual", "description")
View(eset_df) #All Column IDs restored to original 

#2: Filter rows with only "responder" or "non-responder"phenotypes into a new data frame
filtered_df <- eset_df[grepl("responder|nonresponder", eset_df$description) & 
                         !grepl("not assessable", eset_df$description), ]
# This data frame has removed 12 rows or samples which had a phenotype of "non assessable"
View(filtered_df)

# Confirm if there is "not assessable" phenotype at the end of the filtered data 
not_assessable_at_end <- grepl("not assessable$", filtered_df$description)
# Get total sum of the above logical vector
sum_not_assessable_at_end <- sum(not_assessable_at_end)
sum_not_assessable_at_end #This should return 0

############### Data frame cleanup and assigning "responder" or "nonresponder" phenotype annotations to samples in rows #############################

#Remove unwanted columns with non-numeric data except for "description" before differential expression 
# Create a new data frame 
new_df <- filtered_df[, !(names(filtered_df) %in% c("sample", "agent", "individual"))]
View(new_df)

#Annotate phenotype information for each row using description column
#Create a function to detect the presence of "responder" or "nonresponder" 
detect_phenotype <- function(description) {
  if (grepl("responder", description, ignore.case = TRUE) && !grepl("nonresponder", description, ignore.case = TRUE)) {
    return("responder")
  } else if (grepl("nonresponder", description, ignore.case = TRUE)) {
    return("nonresponder")
  } else {
    return(NA) # In case neither is found
  }
}

# Apply the function to each element of the description column
phenotypes <- sapply(new_df[, "description"], detect_phenotype)

# Add a new column with detected phneotype information to the data frame
new_df["phenotype"] = phenotypes
View(new_df)
#The phenotype information in description column should match the "responder" or "nonresponder" status in the new phenotype column upon visual inspection

# Create a new data frame with description column deleted
new_df_1 <- new_df[ ,!grepl("description", colnames(new_df))]
View(new_df_1) #This should have phenotype as the last column and description column deleted

#####Differential expression using t -test with equal variance and multiple hypothesis test with Bonferroni and Benjamini Hochberg FDR methods #########

#install stats package if required 
#install.packages("stats")

# Load stats library if required
#library(stats)

# Check for null values 
if (any(is.na(new_df_1))) {
  # Replace null values with zeros
  new_df_1[is.na(new_df_1)] <- 0
}

# Split the data into "responder" and "nonresponder" groups to get an idea about the sample sizes across each phenotypic group; Not used for subsequent t-test
responder_data <- new_df_1[new_df_1$phenotype == "responder", ]
View(responder_data)
print(nrow(responder_data)) #This should yield "responder" sample size
nonresponder_data <- new_df_1[new_df_1$phenotype == "nonresponder", ]
View(nonresponder_data)
print(nrow(nonresponder_data))#This should yield "nonresponder" sample size

#The sample sizes between groups are unequal
#Welch t-test maybe more robust in this condition 
#However, perform subsequent t-tests assuming equal variance as instructed

#Perform a t-test comparing columns(probes) between two groups defined by the phenotype variable
t_test_results <- apply(new_df_1[, -ncol(new_df_1)], 2, function(probe) {
  t_test <- t.test(probe ~ new_df_1$phenotype, var.equal = TRUE)
  return(c(p_value = t_test$p.value, t_stat = t_test$statistic))
})
View(t_test_results)

# Transpose t_test_results and convert into a data frame
t_test_results_df <- as.data.frame(t(t_test_results))

# Compute Bonferroni-corrected p-values
t_test_results_df$Bonferroni_adjusted_p_value <- p.adjust(t_test_results_df$p_value, method = "bonferroni")

# Compute Benjamini-Hochberg FDR values
t_test_results_df$BH_FDR_adjusted_p_value <- p.adjust(t_test_results_df$p_value, method = "BH") #using fdr here should yield the same outcomes; test if required

#Final results of differential expression and p-value adjustments
View(t_test_results_df)

############### Filtering Bonferroni adjusted p values #############################

# Filter rows based on Bonferroni-adjusted p-value with alpha = 0.01
filtered_df_bonferroni <- t_test_results_df[t_test_results_df$Bonferroni_adjusted_p_value < 0.01, ]
View(filtered_df_bonferroni)

# Create a new data frame without BH values
final_bonferroni_df <- data.frame(
  Probe = rownames(filtered_df_bonferroni),
  T_Stat = filtered_df_bonferroni$t_stat.t,
  p_value = filtered_df_bonferroni$p_value,
  Bonferroni_adjusted_p_value = filtered_df_bonferroni$Bonferroni_adjusted_p_value
)

# Display the new dataframe
View(final_bonferroni_df)

# Display the number of significant genes
cat("Number of genes passing Bonferroni correction with alpha = 0.01", ":", nrow(final_bonferroni_df), "\n")

#Write data frame into comma separated value file
write.csv(final_bonferroni_df, file = "DE_t_test_Bonferroni_alpha1%_Easwaran_Meena.csv", row.names = FALSE)

############### Filtering BH adjusted p values #############################

# Filter rows based on Benjamini-Hochberg FDR-adjusted p-value with alpha = 0.01
filtered_df_BH <- t_test_results_df[t_test_results_df$BH_FDR_adjusted_p_value < 0.01, ]
View(filtered_df_BH)

# Create a new data frame without Bonferroni values
final_BH_df <- data.frame(
  Probe = rownames(filtered_df_BH),
  T_Stat = filtered_df_BH$t_stat.t,
  p_value = filtered_df_BH$p_value,
  BH_adjusted_p_value = filtered_df_BH$BH_FDR_adjusted_p_value
)

# Display the new dataframe
View(final_BH_df)

# Display the number of significant genes
cat("Number of genes passing Benjamini-Hochberg FDRcorrection with alpha = 0.01", ":", nrow(final_BH_df), "\n")

#Write data frame into comma separated value file
write.csv(final_BH_df, file = "DE_t_test_Benjamini_Hochberg_FDR_Easwaran_Meena.csv", row.names = FALSE)
