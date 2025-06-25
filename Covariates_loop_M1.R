library(readr)  
library(readxl)
library(writexl)
library(dplyr)
library(data.table)

(colnames(read.delim("Z:/analyses/ATH_EWAS_phtalates/data/CHAMACOS/results/phenols/model1/PACE_phthalates_CHAMACOS_descriptives_29082023.txt", sep="\t")))

setwd("Z:/analyses/ATH_EWAS_phtalates/data")
setwd("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/data")    #for mobaxterm

cohorts <- c("NELA", "MOCEH", "MMIP", "LINA", "INMA", "HOME", "HELIX", "GenR_Next_NYU", "GenR_Next_NIPH", "GENR", "ELEMENT", "CHAMACOS", "BIS") #13 cohorts
#in phenols, only 12
model_folders <- c("model1", "model2", "model3")
column_names <- c("var", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "SD", "freq")
covariate_summary_df <- data.frame(matrix(ncol=9, nrow = 0))
colnames(covariate_summary_df) <- column_names


for (cohort in cohorts){
  
  for (model in model_folders){
    #select 9y model in case of GENR  
    if (cohort == "GENR" & model == "model2"){
      descriptives_folder <- file.path( ".", cohort, "results", "phenols", model, "mod2_9y")
    }
    else {
    descriptives_folder <- file.path( ".", cohort, "results", "phenols", model)
    }
    
    descriptives_files <- list.files(descriptives_folder, pattern = "descriptives", full.names = TRUE, ignore.case = TRUE)
  
    # Debugging: Print detected files
    print(paste("Files found for", cohort, model, ":", paste(descriptives_files, collapse = ", ")))
    
    
    #filter out INMA BPA descriptives 
    if (cohort == "INMA"){
      descriptives_files <- descriptives_files[!grepl("_bpa_", descriptives_files, ignore.case = TRUE)] 
    }
    #filter out ELEMENT n= 165 descriptives 
    if (cohort == "ELEMENT"){
      descriptives_files <- descriptives_files[!grepl("165", descriptives_files, ignore.case = TRUE)] 
    }
    if (length(descriptives_files) == 0) {
      warning(paste("No descriptives file found for", cohort, model))
      next  # Skip to the next cohort if no file is found
    }
    
    for (descriptives_file in descriptives_files){
      #restructure NELA 
      if (cohort == "NELA"){
        nela_file <- read.csv(descriptives_file, header = TRUE, stringsAsFactors = FALSE)
        cohort_descriptives <- data.frame(matrix(nrow = 28, ncol = 9))
        colnames(cohort_descriptives) <- column_names
        for (i in 1:nrow(nela_file)){
          cleaned_string <- gsub('""', '"', nela_file[i,1])
          split_values <- strsplit(cleaned_string, ",")[[1]]
          rowname <- split_values[1]
          cohort_descriptives[i, ] <- split_values[2:10]
          rownames(cohort_descriptives)[i] <- rowname
        }
        
      } else if (cohort == "HELIX"){
        cohort_descriptives <- read.csv(descriptives_file, header = TRUE, stringsAsFactors = FALSE)
        cohort_descriptives$freq <- as.character(cohort_descriptives$freq)
        print(colnames(cohort_descriptives))
        colnames(cohort_descriptives)[colnames(cohort_descriptives) == "X1st.Qu."] <- "1st Qu."
        colnames(cohort_descriptives)[colnames(cohort_descriptives) == "X3rd.Qu."] <- "3rd Qu."
        }
        
       else {
        cohort_descriptives <- tryCatch({
        read_delim(descriptives_file, delim = "\t", col_types = cols())  # Try tab-separated first
      }, error = function(e) {
        tryCatch({
          fread(descriptives_file, sep =  "","" )  # Try comma-separated if tab fails
        }, error = function(e) {
          warning(paste("Could not read file:", descriptives_file))
          return(NULL)
        })
      })
      }
      
      print(str(cohort_descriptives))
      #convert 
      cols_to_convert <- c("Min.", "Median", "Mean", "1st Qu.", "3rd Qu.", "Max.", "SD")
      
      cohort_descriptives[cols_to_convert] <- lapply(cohort_descriptives[cols_to_convert], function(x) gsub('""', '"', x))
      cohort_descriptives[cols_to_convert] <- lapply(cohort_descriptives[cols_to_convert], function(x) gsub('"', '', x))
      cohort_descriptives[cols_to_convert] <- lapply(cohort_descriptives[cols_to_convert], function(x) gsub(' ', '', x))
      
      # Convert to numeric and handle non-numeric values (like "NA")
      covariate_summary_df[cols_to_convert] <- lapply(covariate_summary_df[cols_to_convert], function(x) as.numeric(as.character(x)))
      cohort_descriptives[cols_to_convert] <- lapply(cohort_descriptives[cols_to_convert], function(x) as.numeric(as.character(x)))
      
      # Convert character columns to character in both data frames
      covariate_summary_df$var <- as.character(covariate_summary_df$var)
      covariate_summary_df$freq <- as.character(covariate_summary_df$freq)
      cohort_descriptives$var <- as.character(cohort_descriptives$var)
      cohort_descriptives$freq <- as.character(cohort_descriptives$freq)

      #make new column with cohort name in it
      cohort_descriptives$cohort <- paste(cohort, model, sep = "_")
      
      if (cohort == "HELIX"){
        if (length(grep(pattern = "HELIX_Descriptives1_", descriptives_file)) > 0) {
          cohort_descriptives$cohort <- paste(cohort, model, "REST", sep= "_")}
        if (length(grep(pattern = "HELIX_Descriptives2_", descriptives_file)) > 0){
          cohort_descriptives$cohort <- paste(cohort, model, "MEPA & ETPA", sep= "_")}
        if (length(grep(pattern = "HELIX_Descriptives3_", descriptives_file)) > 0){
          cohort_descriptives$cohort <- paste(cohort, model, "OXBE", sep= "_")}
      }
      print(paste("Appending", nrow(cohort_descriptives), "rows from", cohort, model, "to final dataframe"))
      
      covariate_summary_df <- bind_rows(covariate_summary_df, cohort_descriptives)
    }
  }
} 
#remove variables related to batch 
rows_to_remove <- grep("batch|plate", covariate_summary_df$var, ignore.case = TRUE)
length(rows_to_remove)
# Remove those rows
covariate_summary_df <- covariate_summary_df[-rows_to_remove, ]
str(covariate_summary_df)
View(covariate_summary_df)

#now split into M1 M2 M3

covariate_summary_df_sorted <- covariate_summary_df[order(-grepl("model1", covariate_summary_df$cohort), 
                         -grepl("model2", covariate_summary_df$cohort),
                         -grepl("model3", covariate_summary_df$cohort),
                         covariate_summary_df$cohort),]
View(covariate_summary_df_sorted)

#now filter based on like terms 
covariate_summary_final <- covariate_summary_df_sorted[order(grepl("edadm", covariate_summary_df_sorted$var, ignore.case = TRUE),  
                                                    grepl("imcm", covariate_summary_df_sorted$var, ignore.case = TRUE),
                                                    grepl("sges", covariate_summary_df_sorted$var, ignore.case = TRUE),
                                                    grepl("NK", covariate_summary_df_sorted$var, ignore.case = TRUE),
                                                    grepl("BCell", covariate_summary_df_sorted$var, ignore.case = TRUE),
                                                    grepl("CD4T", covariate_summary_df_sorted$var, ignore.case = TRUE), 
                                                    grepl("Gran", covariate_summary_df_sorted$var, ignore.case = TRUE),
                                                    grepl("CD8T", covariate_summary_df_sorted$var, ignore.case = TRUE),
                                                    grepl("Mono", covariate_summary_df_sorted$var, ignore.case = TRUE),
                                                    grepl("nRBC", covariate_summary_df_sorted$var, ignore.case = TRUE),
                                                    grepl("parity", covariate_summary_df_sorted$var, ignore.case = TRUE),
                                                    grepl("sexo", covariate_summary_df_sorted$var, ignore.case = TRUE),
                                                    grepl("estudios", covariate_summary_df_sorted$var, ignore.case = TRUE),
                                                    grepl("smok", covariate_summary_df_sorted$var, ignore.case = TRUE),
                                                    grepl("BMI", covariate_summary_df_sorted$var, ignore.case = TRUE),
                                                    covariate_summary_df_sorted$var),]



write_xlsx(covariate_summary_final, "Phenols_Covariate_Summary.xlsx")
  
