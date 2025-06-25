
# =========================================================================== #
#sylviefraley
#Model3, cross sectional childhood
# =========================================================================== #

# =========================================================================== #
# The following R code will allow you to complete the EWAS Prenatal phtalates
# and bisphenols VS child blood included in the PACE phtalates/bisphenols and
# child blood DNA methylation analysis plan.
#
# The code also produces files summarizing the variables included in the EWAS.
# You shouldn't have to rewrite or add anything to the following code, unless
# otherwise stated.
#
# There are just two inputs required for this analysis:
# 1) pheno: a dataframe containing all the "phenotype" data needed for this
#    project. Each row is a sample (individual) and each column is a different
#    variable. 
#    - Exposure:
#       --- Phtalates ---
#       * "MEP_mraw_creat_preg_log2"
#       * "MIBP_mraw_creat_preg_log2"
#       * "MNBP_mraw_creat_preg_log2"
#       * "MBZP_mraw_creat_preg_log2"
#       * "MEHP_mraw_creat_preg_log2"
#       * "MEHHP_mraw_creat_preg_log2"
#       * "MEOHP_mraw_creat_preg_log2"
#       * "MECPP_mraw_creat_preg_log2"
#       ---Bisphenols ---
#       * "BPA_mraw_creat_preg_log2"
#       * "BPS-BPF_mraw_creat_preg_log2"
#       * "PARABENS_mraw_creat_preg_log2"
#       * "TRICLOSAN_mraw_creat_preg_log2"
#       * "BP_mraw_creat_preg_log2"
#    - Main covariates:
#       * "SampleID"
#       * "sexo" (Sex)
#       * "edadm" (Maternal age)
#       * "estudios3c" (Maternal education)
#       * "smok_preg" (3-level maternal pregnancy smoking status)
#       * "zBMI" (child BMI as z score)
#       * "age" (child's age)
#       * "gestage" (Gestational Age)
#       * "NK"
#       * "Gran",
#       * "Bcell"
#       * "CD8T"
#       * "CD4T"
#       * "Mono"
#    - Optional covariates:
#       * "anc" (ancestry)
#       * "batch" (technical covariates)
#       * "sel" (selection factor)
#
# If these columns are named differently in your dataset, please rename the
# columns accordingly. Details on how to code these variables are provided
# in the analysis plan.
#
# 2) beta_matrix: a matrix of methylation illumina beta values. Each column is
#    a sample and each row is a probe on the array (450k or EPIC). 
#    Column names must correspond to the SampleID in pheno.
#
# =========================================================================== #

# =========================================================================== #
### Go to the directory where to save result
### (change to your own working directory)
#setwd("Z:/analyses/ATH_EWAS_phenols/results")
#setwd("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phenols/results") #mobaxterm

### Install packages if they are not already installed
install.packages("haven")
if (!require("openxlsx", quietly = TRUE))
  install.packages("openxlsx")
if (!require("readxl", quietly = TRUE))
  install.packages("readxl")
if (!require("data.table", quietly = TRUE))
  install.packages("data.table")
if (!require("sandwich", quietly = TRUE))
  install.packages("sandwich")
if (!require("lmtest", quietly = TRUE))
  install.packages("lmtest")
if (!require("parallel", quietly = TRUE))
  install.packages("parallel")
if (!require("R.utils", quietly = TRUE))
  install.packages("R.utils")
if (!require("Hmisc", quietly = TRUE))
  install.packages("Hmisc")
if (!require("psych", quietly = TRUE))
  install.packages("psych")
if (!require("minfi", quietly = TRUE))
  install.packages("minfi")
if (!require("xlsx", quietly = TRUE))
  install.packages("xlsx")
if (!require("qqman", quietly = TRUE))
  install.packages("qqman")
if (!require("matrixStats", quietly = TRUE))
  install.packages("matrixStats")
if (!require("writexl", quietly = TRUE))
  install.packages("writexl")
### Load required packages 
library(haven)
library(dplyr)
library(openxlsx) # to open excel
library(readxl) # to get data out of excel and into R
library(data.table) # to process results
library(sandwich) # to estimate the standard error
library(lmtest) # to use coeftest
library(parallel) # to use multicore approach - part of base R
library(R.utils) # utility functions
library(Hmisc) # to describe function
library(psych) # to run correlation analysis
library(minfi) # to get DNA methylation data
library(xlsx) # to open excel
library(qqman) # to get QQ plots
library(matrixStats)
library(writexl)

# Load and check phenotype data ============================================= #

# pheno: a dataframe containing all the "phenotype" data needed for this project. 
# Each row is a sample (individual) and each column is a different variable.
# Ensure all traits and covariates have been derived as specified in the
# analysis plan and that exclusions have been made.

# Remember not to include any ehtnicities with less than 20 cases in the cohort.
# If your study has several ethnicities at a substantial frequency, then perform
# an analysis with ALL the participants adjusting for the best ethnic variable
# available (read protocol for more details)

# Read pheno file (change to your file name and directory)
# pheno <- foreign::read.spss("/data/RAW_DATA/450K/SAB_450K_blood/EWAS_Phtalates_Metil_cord_Sab.sav",
#                            to.data.frame = TRUE)

#do imputation and zBMI 
#setwd("Z:/analyses/ATH_EWAS_phenols/db/")
setwd("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phenols/db") #mobaxterm

pheno <- read.table('../db/pheno_child/pheno_ALL.txt', sep = '\t', header = T, stringsAsFactors = T)
dim(pheno)#1114


raw_values <- read_dta("./biomarker_final_database.dta")
dim(raw_values) #1330 IDs, 371 columns
names(raw_values)


#remove NAs for creatinine from raw_values
no_NA_creatinine <- raw_values[!is.na(raw_values$hs_creatinine_cg), ] #29 NAs removed
dim(no_NA_creatinine) #1301
colnames(no_NA_creatinine)

#LODS/2
lods <- c(
  bpa = 0.03 / 2,
  mepa = 0.03 / 2
)

#looking for descriptives of compounds of interest in CHILDREN 
id_cols <- c("cohort", "helix_id", "sample_id", "subject_id", "hs_creatinine_cg", "hs_creatinine_cdesc")
matching_cols <- names(no_NA_creatinine)[
  grepl("bpa_c(desc)?$|mepa_c(desc)?$", names(no_NA_creatinine), ignore.case = TRUE)
]
keep_cols <- c(id_cols, matching_cols)
subset_desc <- no_NA_creatinine[, keep_cols]

# 2. Filter rows where any "desc" column 
desc_cols <- grep("desc", names(subset_desc), ignore.case = TRUE, value = TRUE)
print(desc_cols)

descriptives <- data.frame(
  chemical = c("bpa", "mepa"),
  N_below_LOD = 0,
  min = NA, IQR = NA, p5 = NA, p25 = NA, p50 = NA, p75 = NA, p95 = NA,
  max = NA, mean = NA, SD = NA,
  stringsAsFactors = FALSE
)
rownames(descriptives) <- descriptives$chemical
print(descriptives)

#replace desc = 2 NAs with LOD, no other changes
count_imputed <- 0
for (i in 1:nrow(subset_desc)){
  for(chem in names(lods)){
    level_col <- paste0("hs_", chem, "_c")
    desc_col  <- paste0("hs_", chem, "_cdesc")
    
    # Check if desc = 2, if so change imputaiton 
    if (subset_desc[i, desc_col] == (2)) {
      subset_desc[i, level_col] <- (lods[chem])
      count_imputed <- count_imputed + 1
      if (chem %in% rownames(descriptives)){
        descriptives[chem, "N_below_LOD"] <- descriptives[chem, "N_below_LOD"] + 1
      }
    }
  }
  
}
print(count_imputed)
print(descriptives)


#now, normalize by creatinine and log2 transform for all values in desc 1 and 2 
BPA_normed <- numeric()
MEPA_normed <- numeric()
for (i in 1:nrow(subset_desc)){
  for (chem in names(lods)) {
    level_col <- paste0("hs_", chem, "_c")
    desc_col <- paste0("hs_", chem, "_cdesc")
  
    if (level_col %in% names(subset_desc) && !is.na(subset_desc[i,level_col]) && (subset_desc[i, desc_col]) %in% c(1, 2)) {
    normed <- as.numeric(subset_desc[i, level_col]) / as.numeric(subset_desc[i, "hs_creatinine_cg"])
      if (subset_desc[i, desc_col] == 1){
        if (chem == "bpa"){
          BPA_normed <- c(BPA_normed, normed)
        }
        if (chem == "mepa"){
          MEPA_normed <- c(MEPA_normed, normed)
        }
      }
    if (!is.na(normed) && normed > 0) {
      subset_desc[i, level_col] <- log2(normed)
    } else {
      print("issue with creatine normalization")
      subset_desc[i, level_col] <- NA  # or leave it as-is
    }
    }
  }
}
summary(subset_desc$hs_bpa_c) #now negative values bc log2 transformed

#create descriptives with values not log2 transformed
descriptives["bpa", "min"] <- min(BPA_normed, na.rm = TRUE)
descriptives["mepa", "min"] <- min(MEPA_normed, na.rm = TRUE)

descriptives["bpa", "p5"] <- quantile(BPA_normed, 0.05, na.rm = TRUE)
descriptives["mepa", "p5"] <- quantile(MEPA_normed, 0.05, na.rm = TRUE)

descriptives["bpa", "p25"] <- quantile(BPA_normed, 0.25, na.rm = TRUE)
descriptives["mepa", "p25"] <- quantile(MEPA_normed, 0.25, na.rm = TRUE)

descriptives["bpa", "p50"] <- quantile(BPA_normed, 0.5, na.rm = TRUE)
descriptives["mepa", "p50"] <- quantile(MEPA_normed, 0.5, na.rm = TRUE)

descriptives["bpa", "p75"] <- quantile(BPA_normed, 0.75, na.rm = TRUE)
descriptives["mepa", "p75"] <- quantile(MEPA_normed, 0.75, na.rm = TRUE)

descriptives["bpa", "p95"] <- quantile(BPA_normed, 0.95, na.rm = TRUE)
descriptives["mepa", "p95"] <- quantile(MEPA_normed, 0.95, na.rm = TRUE)

descriptives["bpa", "max"] <- max(BPA_normed, na.rm = TRUE)
descriptives["mepa", "max"] <- max(MEPA_normed, na.rm = TRUE)

descriptives["bpa", "mean"] <- mean(BPA_normed, na.rm = TRUE)
descriptives["mepa", "mean"] <- mean(MEPA_normed, na.rm = TRUE)

descriptives["bpa", "SD"] <- sd(BPA_normed, na.rm = TRUE)
descriptives["mepa", "SD"] <- sd(MEPA_normed, na.rm = TRUE)

descriptives["bpa", "IQR"] <- IQR(BPA_normed, na.rm = TRUE)
descriptives["mepa", "IQR"] <- IQR(MEPA_normed, na.rm = TRUE)

print(descriptives)
write_xlsx(descriptives, file.path(getwd(), "descriptives_M3.xlsx"))

#now correlation 
# At this point, subset_desc has already done imputation and creatinine transformation - again making it for real values so desc 1 and 2. 
bpa_values <- c()
mepa_values <- c()
for (i in 1:nrow(subset_desc)){
  if (subset_desc[i, "hs_bpa_cdesc"] %in%  c(1,2)  && subset_desc[i, "hs_mepa_cdesc"] %in%  c(1,2) ) {
    bpa_values <- c(bpa_values, subset_desc[i, "hs_bpa_c"])
    mepa_values <- c(mepa_values,subset_desc[i, "hs_mepa_c"])
  }

}
correlation <- matrix(nrow = 2, ncol = 2)
rownames(correlation) <- c("BPA", "MEPA")
colnames(correlation) <- c("BPA", "MEPA")
correlation[1,2] <- cor(as.numeric(bpa_values), as.numeric(mepa_values), use = "complete.obs")
correlation[2, 1] <- cor(as.numeric(mepa_values), as.numeric(bpa_values), use = "complete.obs")
correlation[1, 1] <- cor(as.numeric(bpa_values), as.numeric(bpa_values), use = "complete.obs")
correlation[2, 2] <- cor(as.numeric(mepa_values), as.numeric(mepa_values), use = "complete.obs")

getwd()
correlation <- data.frame(correlation)
write_xlsx(correlation, file.path(getwd(), "correlations_M3.xlsx"))

##-----------------------------------------------------------------
#compare mine to pheno 
  
names(pheno)[names(pheno) == "HelixID"] <- "helix_id"

pheno_col <- c(
  "BPA_mraw_creat_preg_log2", "PARABENS1_MEPA_mraw_creat_preg_log2")


common_ids <- intersect(subset_desc$helix_id, pheno$helix_id)
length(common_ids) #1113 IDs common to both dataframes

subset_desc_aligned <- subset_desc[subset_desc$helix_id %in% common_ids, ]
subset_desc_aligned <- subset_desc_aligned[order(subset_desc_aligned$helix_id), ]
pheno_aligned <- pheno[pheno$helix_id %in% common_ids, ]
pheno_aligned <- pheno_aligned[order(pheno_aligned$helix_id), ]

comparison_list <- list()

for (chem in names(lods)) {
  level_col <- paste0("hs_", chem, "_c")
  pheno_match <- grepl(chem, pheno_col, ignore.case = TRUE)
  pheno_col_this <- pheno_col[pheno_match]
  
  if (length(pheno_col_this) == 1 && level_col %in% names(subset_desc_aligned) && pheno_col_this %in% names(pheno_aligned)) {
    merged_df <- merge(
      subset_desc_aligned[, c("helix_id", level_col)],
      pheno_aligned[, c("helix_id", pheno_col_this)],
      by = "helix_id"
    )
    
    # Make sure there's data to summarize
    if (nrow(merged_df) > 0) {
      x <- merged_df[[level_col]]
      y <- merged_df[[pheno_col_this]]
      
      if (all(is.na(x)) || all(is.na(y))) {
        corr <- NA
        sum_x <- rep(NA, 6)
        names(sum_x) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
        sum_y <- sum_x
      } else {
        sum_x <- summary(x)
        sum_y <- summary(y)
        corr <- cor(x, y, use = "complete.obs")
      }
      
      comparison_row <- data.frame(
        Chemical = chem,
        NA_Min = as.numeric(sum_x["Min."]),
        NA_Q1 = as.numeric(sum_x["1st Qu."]),
        NA_Median = as.numeric(sum_x["Median"]),
        NA_Q3 = as.numeric(sum_x["3rd Qu."]),
        NA_Max = as.numeric(sum_x["Max."]),
        Pheno_Min = as.numeric(sum_y["Min."]),
        Pheno_Q1 = as.numeric(sum_y["1st Qu."]),
        Pheno_Median = as.numeric(sum_y["Median"]),
        Pheno_Q3 = as.numeric(sum_y["3rd Qu."]),
        Pheno_Max = as.numeric(sum_y["Max."]),
        Correlation = corr,
        stringsAsFactors = FALSE
      )
      
      comparison_list[[length(comparison_list) + 1]] <- comparison_row
    }
  }
}

# Only combine non-null rows
comparison_list <- Filter(Negate(is.null), comparison_list)
comparison <- do.call(rbind, comparison_list)
print(comparison)
write_xlsx(comparison, "imputation_comparison_M3.xlsx")

# --- Set initial parameters ------------------------------------------------ #

#import new values into pheno data frame 

names(pheno)[names(pheno) == "BPA_mraw_creat_preg_log2"] <- "hs_bpa_c"
names(pheno)[names(pheno) == "PARABENS1_MEPA_mraw_creat_preg_log2"] <- "hs_mepa_c"
update_cols <- c("hs_bpa_c", "hs_mepa_c")
id_col_1 <- "helix_id"

update_rows_by_id <- function(pheno, subset_desc, id_col_1, update_cols) {
  # Ensure ID column is character for safe merging
  pheno[[id_col_1]] <- as.character(pheno[[id_col_1]])
  subset_desc[[id_col_1]] <- as.character(subset_desc[[id_col_1]])
  
  # Check update_cols exist in both dataframes
  missing_cols <- setdiff(update_cols, intersect(names(pheno), names(subset_desc)))
  if (length(missing_cols) > 0) {
    stop("These update_cols are missing in one or both dataframes: ", paste(missing_cols, collapse = ", "))
  }
  
  # Match rows based on ID
  matched_indices <- match(pheno[[id_col_1]], subset_desc[[id_col_1]])
  
  # Update values only where ID matches
  for (col in update_cols) {
    pheno[[col]][!is.na(matched_indices)] <- subset_desc[[col]][matched_indices[!is.na(matched_indices)]]
  }
  
  return(pheno)
}
update_rows_by_id(pheno, subset_desc, id_col_1, update_cols)


cohort <- "HELIX" # Define cohort name
samplename <- "helix_id" # Define Sample ID column

# Define exposure traits
# This is important to define traits so that the loop to execute EWAS works.

traits <- c("hs_bpa_c", "hs_mepa_c")


# Define CATEGORICAL VARIABLES

cat_vars <- c('parity', 'sexo', 'estudios3c', 'smok_preg')
length(cat_vars)

# Define CONTINUOUS VARIABLES, changed this from prepreg to z
cont_vars <- c('edadm', 'zBMI', 'age', 'gestage', 'NK', 'Bcell', 'CD4T', 'Gran', 'CD8T', 'Mono',
               'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
length(cont_vars)

# Covariates "batch", "and" and "sel" are optional, you can include them
#did not see batch as an option, skipped this
# in the corresponding vector (categorical / continuous)
#batch_ones <- (grepl( "batch", names(raw_values)))
#subset_batch <- raw_values[,batch_ones]
#View(subset_batch)
# Preferred categorization of smok_preg is into 3 groups

# Date in day/month/year for filenames
datefn <- format(Sys.Date(), "%d%m%Y")

# Prefix path for filenames
prefix <- './PACE_phenols_newimp_'

# Create directory for the results
dir.create('M3_cross_sectional_chilhood_newimputation')


# Check if all needed variables are present in phenotype file --------------- #
# If not, please correct.
phenocols <- c(samplename, traits, cat_vars, cont_vars)
print(phenocols)
for (col in phenocols) {
  if (col %in% colnames(pheno)) {
    print(paste("variable called", col, "is present in pheno"))
  } else {
    print(paste("CAUTION: the variable", col, "is missing from pheno"))
  }
}


# Row names as sample ID
rownames(pheno) <- pheno$SampleID

# Make sure variables are defined properly ---------------------------------- #

# Converting continuous variables to numeric class
pheno[, cont_vars] <- sapply(pheno[, cont_vars], as.numeric)
sapply(pheno[, cont_vars], class) # Check

# Converting categoricals variables to factor class
pheno[, cat_vars] <- lapply(pheno[, cat_vars], factor)
sapply(pheno[, cat_vars], class) # Check

# Round decimals
pheno$zBMI <- round(pheno$zBMI, digits = 3)
pheno$gestage <- round(pheno$gestage, digits = 1)
pheno$NK <- round(pheno$NK, digits = 4)
pheno$Bcell <- round(pheno$Bcell, digits = 4)
pheno$CD4T <- round(pheno$CD4T, digits = 4)
pheno$Gran <- round(pheno$Gran, digits = 4)
pheno$CD8T <- round(pheno$CD8T, digits = 4)
pheno$Mono <- round(pheno$Mono, digits = 4)

# Subsetting dataset to the selected variables
phenocols %in% colnames(pheno)

pheno <- pheno[,phenocols] #still 1114 rows but now subsetting columns 
dim(pheno) #27 columns, same number of rows (1114)

# check missings in phenols

summary(pheno$hs_mepa_c, useNA="ifany")#2 NA
summary(pheno$hs_bpa_c, useNA="ifany")#11 NA

# Select complete cases to do the analysis, check this matches descriptives
# If you have >5% missing for any of the variables, check with us before proceeding
pheno <- pheno[complete.cases(pheno), ] #remove rows with any NA values in any column
dim(pheno) #now 1024 rows now, 
colnames(pheno)

# Load methylation data ===================================================== #

# beta_matrix: a matrix of methylation illumina beta values. Each column is a
# sample and each row is a probe on the array (450k or EPIC). 
# Column names must correspond to the Sample_ID in pheno.

# Load methylation data
#..# load('/data/RAW_DATA/450K/SAB_450K_blood/SAB_Blood_0y_450K.Rdata')
load("/PROJECTES/HELIX_OMICS/data_final/methyl/child/8y/blood_450K_QChelix_20170401/methylome_subcohort_ComBatSlide_6cells_notfitr_v4.Rdata")

# Extract beta matrix from the GenomicRatioSet
beta_matrix <- beta_matrix2 <- getBeta(methylome_subcohort_ComBatSlide_notfitr)
dim(beta_matrix) #480,444 rows, 1192 columns

# =========================================================================== #
### Select subset of children for the analysis and match IDs

##########################################

# bpa,mepa

# Filter beta_matrix to have the same individuals as in the pheno file
beta_matrix <- beta_matrix2[, colnames(beta_matrix2) %in% rownames(pheno)]
dim(beta_matrix) #now 1024 columns to match the samples in pheno1

common_ids <- intersect(rownames(pheno), colnames(beta_matrix))
pheno <- pheno[common_ids, ]
beta_matrix <- beta_matrix[, common_ids]

# Match
beta_matrix <- beta_matrix[, match(rownames(pheno), colnames(beta_matrix))]

# Match IDs in beta_matrix with pheno file
table(ifelse(rownames(pheno) == colnames(beta_matrix),
             "MATCHED", "Not Matched"))
#confirmed all matched

# Save new pheno file (change output file name)
write.csv(pheno,
          paste0(prefix, cohort, '_phenoNEWIMP_', datefn, '.csv'),
          row.names = F)


# =========================================================================== #
### Descriptive of pheno data

n <- nrow(pheno)

summarize_categ <- function(categorical) {
  t_categorical <- table(pheno[[categorical]])
  names(t_categorical) <- paste(categorical, names(t_categorical))
  return(as.data.frame(t_categorical))
}

cat_descriptive <- lapply(cat_vars, summarize_categ)
cat_descriptive <- do.call(rbind, cat_descriptive)
colnames(cat_descriptive) <- c('var', 'freq')

sums <- apply(pheno[, cont_vars], 2, summary) # Summaries
sds <- apply(pheno[, cont_vars], 2, sd) # Standard deviations
cont_descriptive <- t(rbind(sums, 'SD' = sds))
cont_descriptive <- round(cont_descriptive, 2)
cont_descriptive <- as.data.frame(cbind('var' = rownames(cont_descriptive),
                                        cont_descriptive))
rownames(cont_descriptive) <- NULL

table_descriptive <- merge(cont_descriptive, cat_descriptive, all = T, sort = F)


# Save results
write.csv(table_descriptive,
          paste0(prefix, cohort, "_DescriptivesNEWIMP_", datefn, ".csv"),
          row.names = F)


# =========================================================================== #
# Winsorization method to remove outliers (check if this has not already been
# applied to your data)

# Check that CpGs are rows
beta_matrix[1:5, 1:10] #yes

# Function
winsorize <- function(methylation, pct) {
  
  stopifnot(is.matrix(methylation))
  if (nrow(methylation) < ncol(methylation)) {
    warning("expecting CpG methylation as rows (long dataset)")
  }
  
  quantiles <- matrixStats::rowQuantiles(methylation,
                                         probs = c(pct, 1 - pct),
                                         na.rm = T)
  low <- quantiles[, 1]
  upper <- quantiles[, 2]
  
  outliers.lower <- rowSums(methylation < low, na.rm = T)
  outliers.upper <- rowSums(methylation > upper, na.rm = T)
  
  idx <- which(methylation < low, arr.ind = T)
  methylation[idx] <- low[idx[, 1]]
  
  idx <- which(methylation > upper, arr.ind = T)
  methylation[idx] <- upper[idx[, 1]]
  
  n <- rowSums(!is.na(methylation))
  log <- data.frame(outliers.lower, outliers.upper, n)
  
  return(list(methylation = methylation, log = log))
}

# Replace outliers
replace_outliers <- winsorize(beta_matrix, 0.005)
beta_matrix <- as.matrix(replace_outliers$methylation)

# Save winsorization log
wins_log <- as.data.frame(replace_outliers$log)

write.table(wins_log, file = paste0(prefix, cohort, "_PACE_INTEXT_winslog_NEWIMP_", datefn, ".txt"),
            sep = "\t", col.names = T, row.names = T, append = F, quote = F)


# ========================================================================== #
### Initial checks
# Check that rows are samples and columns are probes on the array (450k or
# EPIC). If not, then transpose betas so that rows are samples and columns are
# probes on the array (450k or EPIC)
dim(beta_matrix)
beta_matrix <- t(beta_matrix)
dim(beta_matrix)

# Check that all SampleIDs are the same - need to be in the same order!
stopifnot(rownames(beta_matrix) == rownames(pheno))


# Preparing Functions
apply_model <- function(beta, phen, form) {
  phen$CpGi <- beta
  coeff <- coef(summary(lm(formula = form, data = phen)))
  return(list(estimate <- coeff[, 1],
              SE <- coeff[, 2],
              tval <- coeff[, 3],
              pval <- coeff[, 4],
              N <- length(which(!is.na(phen$CpGi))) )
  )
}

get_model_data <- function(element, data, nrows) {
  mat <- matrix(unlist(lapply(data, "[[", element)), nrow = nrows, byrow = TRUE)
  if (dim(mat)[2] > 1) {
    colnames(mat) <- names(data[[element]][[1]])
  }
  rownames(mat) <- names(data)
  return(mat)
}

# =========================================================================== #
# Main loop for correlation with blood cell proportions, EWAS (different models
# per variable), lambda calculation and QQ-plot

traits <- c("hs_bpa_c",
            'hs_mepa_c')
# bpa,mepa

for (trait in traits) {
  
  # Association between blood cell type proportion and phthalates/bisphenols
  NK <- lm(pheno$NK ~ pheno[[trait]])
  NK <- summary(NK)$coefficients[2, ]
  Bcell <- lm(pheno$Bcell ~ pheno[[trait]])
  Bcell <- summary(Bcell)$coefficients[2, ]
  CD4T <- lm(pheno$CD4T ~ pheno[[trait]])
  CD4T <- summary(CD4T)$coefficients[2, ]
  Gran <- lm(pheno$Gran ~ pheno[[trait]])
  Gran <- summary(Gran)$coefficients[2, ]
  CD8T <- lm(pheno$CD8T ~ pheno[[trait]])
  CD8T <- summary(CD8T)$coefficients[2, ]
  Mono <- lm(pheno$Mono ~ pheno[[trait]])
  Mono <- summary(Mono)$coefficients[2, ]
  
  rb <- rbind(NK, Bcell, CD4T, Gran, CD8T, Mono)
  rb <- as.data.frame(rb)
  rb$celltype <- rownames(rb)
  rb$measure <- trait
  rb <- rb[, c(6, 5, 1:4)]
  rownames(rb) <- seq_along(rb)
  
  write.table(rb,
              file = paste0(prefix, cohort, "_LMCellProp_", trait, '_', datefn, ".tsv"),
              row.names = F, quote = F, sep = '\t')
  
  # Remove unnecessary variables
  rm(list = c('rb', 'NK', 'Bcell', 'CD4T', 'Gran', 'CD8T', 'Mono'))
  gc()
  
  
 
  #  covariate names in your pheno data 
  
  ModelA <- as.formula(paste0('CpGi ~ ', trait, ' + edadm + zBMI + parity + estudios3c + smok_preg + sexo + age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10'))
  ModelB <- as.formula(paste0('CpGi ~ ', trait, ' + edadm + zBMI + parity + estudios3c + smok_preg + sexo + age + Mono + Bcell + CD4T + CD8T +  Gran + NK + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10'))
  ModelC <- as.formula(paste0('CpGi ~ ', trait, ' + edadm + zBMI + parity + estudios3c + smok_preg + sexo + age + Mono + Bcell + CD4T + CD8T +  Gran + NK + gestage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10'))
  

  
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  ncpgs <- ncol(beta_matrix)
  models <- c('ModelA', 'ModelB', 'ModelC')
  #duplicate logic, models as strings instead
  formulas <- c(paste0('CpGi ~ ', trait, ' + edadm + zBMI + parity + estudios3c + smok_preg + sexo + age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10'),
                paste0('CpGi ~ ', trait, ' + edadm + zBMI + parity + estudios3c + smok_preg + sexo + age + Mono + Bcell + CD4T + CD8T +  Gran + NK + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10'),
                paste0('CpGi ~ ', trait, ' + edadm + zBMI + parity + estudios3c + smok_preg + sexo + age + Mono + Bcell + CD4T + CD8T +  Gran + NK + gestage + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10'))
  
  
  
  sapply(seq_along(models), function(i) {
    
    results <- apply(beta_matrix, 2, apply_model, phen = pheno, form = as.formula(formulas[i]))
    
    # Get data from results
    beta <- get_model_data(1, results, ncpgs) ; beta <-  beta[, trait]
    SE <- get_model_data(2, results, ncpgs) ;  SE <-  SE[, trait]
    pval <- get_model_data(4, results, ncpgs); pval <-  pval[, trait]
    N <- get_model_data(5, results, ncpgs); N <-  N[, 1]
    
    rm(results); gc() # Remove unnecessary variables
    
    write.table( data.frame(beta = beta, SE = SE, pval = pval, N = N),
                 paste0(prefix, cohort, '_EWAS_', models[i],  '_', trait, '_', datefn, ".tsv"),
                 quote = F, sep = '\t')
    
    
    rm(list = c( 'beta', 'SE', 'N' )) ; gc()   # Remove unnecessary variables
    
    # Calculate lambdas
    lambda <- qchisq( median( pval, na.rm = T), df = 1, lower.tail = F) / qchisq(0.5, 1)
    
    write.table(lambda,
                paste0(prefix, cohort, "_Lambda_", models[i], '_', trait, '_', datefn, ".tsv"),
                col.names = F, quote = F, sep = '\t')
    
    # QQPlot
    pdf(paste0(prefix, cohort, "_QQPlot_", models[i], '_', trait, '_', datefn, ".pdf"))
    qq( pval, main = paste("PACE_phenols", cohort, models[i], trait, sep = ' '))
    dev.off()
    
    rm(list = c( 'pval', 'lambda')); gc()   # Remove unnecessary variables
    
    
  })
}

#if QQplots have trouble generating
setwd("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phenols/db/M3_cross_sectional_chilhood_newimputation")
pval_files <- list.files(pattern = "EWAS_.*\\.tsv")

for (val in pval_files){
  name <- sub("^.*EWAS_", "", val)
  name <- sub("\\.tsv$", "", name)
  
  # Build output filename and plot title
  outfile <- paste0("PACE_phenols_newimp_HELIX_QQPlot_", name, ".pdf")
  plottitle <- paste("QQ Plot:", name)
  ewas <- read.table(val, header = TRUE)
  
  # Filter pvals
  pval <- ewas$pval
  pval <- pval[!is.na(pval) & pval > 0 & pval <= 1]
  
  # Plot
  pdf(outfile)
  qq(pval, main = plottitle)
  dev.off()
}


