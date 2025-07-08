##############################################################
# phenols EWAS EASIER QC
# sylvie fraley
# 13.06.2025
#model 1 - thesis run without INMA BPA, CHAMACOS, GENRNEXT, MMIP, MOCEH, HOME
# vignette: https://github.com/isglobal-brge/EASIER/blob/main/vignettes/EASIER.html
##############################################################


## ################################################## ##
##  Quality Control Script to use with EASIER package ##
##                                                    ##
##  script version; 0.1.2.27                          ##
## ################################################## ##


## -------------------------------------
##  Install EASIER Package Code
## -------------------------------------
##
##  Uncomment this code to install EASIER package
#
# # Install devtools
# install.packages("devtools")
#
# # Install required packages
# devtools::source_url("https://raw.githubusercontent.com/isglobal-brge/EASIER/HEAD/installer.R")

# # Install EASIER package
# devtools::install_github("isglobal-brge/EASIER@HEAD")

##  END -  Install EASIER Package Code
## -------------------------------------


# load package
library(EASIER)
library(dplyr)
library(Cairo)
library(ggplot2)
library(readr)



########## ----------  VARIABLES DEFINED BY USER  ----------  ##########



# Set working directory to metaanalysis folder
setwd("Z:/analyses/ATH_EWAS_phtalates/data")
setwd("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/data")    #for mobaxterm


files_GenR_0y <- list.files(path = "./GENR/results/phenols/model1", pattern = "PACE_phthalates_GenR_EWAS_Model(A|B|C)_BPAugg_average_imputed_log2_.*", full.names = TRUE)
files_GenR_0y # 3 files

#for now remove BPA
files_INMA <- list.files(path = "./INMA/results/phenols/model1", pattern = "PACE_phthalates_INMA_EWAS_Model(A|B|C)_(trcs|etpa|bp3|mepa|prpa)_mraw_(creat_total_log2)_.*", full.names = TRUE)
files_INMA #15 files, all total, leave out BPA

#######

files_NELA <- list.files(path = "./NELA/results/phenols/model1", pattern = "PACE_phthalates_NELA_EWAS_Model(A|B|C)_(ad.*)_(EPB|MPB|BP3|PPB|BPA)(_media)?_S24_.*", full.names = TRUE)
files_NELA #15 files

files <- c (files_GenR_0y, files_INMA, files_NELA)
files 


# Result folder
results_folder <- 'QC_Results_Phenols_Sfraley_thesis'

# Prefixes for each file - ADD
prefixes <- c(
  'Model1A_GenR_BPA', 'Model1B_GenR_BPA', 'Model1C_GenR_BPA',
  'Model1A_INMA_BP3',  'Model1A_INMA_ETPA', 'Model1A_INMA_MEPA', 'Model1A_INMA_PRPA', 'Model1A_INMA_TCS',
  'Model1B_INMA_BP3',  'Model1B_INMA_ETPA', 'Model1B_INMA_MEPA', 'Model1B_INMA_PRPA', 'Model1B_INMA_TCS',
  'Model1C_INMA_BP3',  'Model1C_INMA_ETPA', 'Model1C_INMA_MEPA', 'Model1C_INMA_PRPA', 'Model1C_INMA_TCS',
  'Model1A_NELA_BP3', 'Model1A_NELA_ETPA', 'Model1A_NELA_MEPA', 'Model1A_NELA_PRPA', 'Model1A_NELA_BPA',
  'Model1B_NELA_BP3', 'Model1B_NELA_ETPA', 'Model1B_NELA_MEPA', 'Model1B_NELA_PRPA', 'Model1B_NELA_BPA',
  'Model1C_NELA_BP3', 'Model1C_NELA_ETPA', 'Model1C_NELA_MEPA', 'Model1C_NELA_PRPA', 'Model1C_NELA_BPA'
)
length(prefixes)


# Exclude - MASK snp5
ethnic<-rep("GMAF1p", times = length(prefixes))

# Array type, used : 
artype <- c('450K', '450K', '450K', #GEnR
            '450K', '450K', '450K', '450K', '450K', '450K', '450K', '450K', '450K','450K', '450K', '450K', '450K', '450K', '450K', #INMA
            'EPIC', 'EPIC', 'EPIC','EPIC', 'EPIC', 'EPIC','EPIC', 'EPIC', 'EPIC','EPIC', 'EPIC', 'EPIC','EPIC', 'EPIC', 'EPIC' #NELA
)
length(artype)

# Parameters to exclude CpGs
exclude <- c( 'control_probes',
              'noncpg_probes',
              'Sex',
              'MASK_mapping',
              'MASK_sub30_copy',
              'MASK_extBase',
              'MASK_typeINextBaseSwitch',
              'MASK_snp5_GMAF1p',
              'Unrel_450_EPIC_blood')

N <- c(
  276, 276, 276, #GenR0y, changed already
  316, 314, 316, 315, 316, 316, 314, 316, 315, 316, 316, 314, 316, 315, 316, #INMA
  427, 427, 427, #LINA, changed already
  314, 314, 314, 314, 314, 314, 314, 314, 314, 314, 314, 314, 314, 314, 314 #NELA
)
length(N)
n <- c(NA)

# Minimum sample representation percentage required for CpGs
# Filter minimum percentage of missings for each CpG in cohort
# We need to define two parameters,
#  - colname_NforProbe:
#        Column name with Number of individuals per probe, this variable only needs
#           to be defined if you want to filter CpGs with low representation.
#         If defined value in colname_NforProbe not exists, no filter will be applied
#  - pcMissingSamples :
#        Máximum percent of missing samples allowed,

colname_NforProbe <- 'N_for_probe'
pcMissingSamples <- 0.9

########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########



## ###################### ##
##  QC - Quality Control  ##
## ###################### ##

# Variable declaration to perform precision plot
medianSE <- numeric(length(files))
value_N <- numeric(length(files))

if(length(n) == length(N))
  value_n <- numeric(length(files))

cohort_label <- character(length(files))

# Prepare output folder for results (create if not exists)
if(!dir.exists(file.path(getwd(), results_folder )))
  suppressWarnings(dir.create(file.path(getwd(), results_folder)))

# Remove tmp files
if( file.exists(paste0(results_folder,"/tmp_pretQC.txt")) )
  file.remove(paste0(results_folder,"/tmp_pretQC.txt"))
if( file.exists(paste0(results_folder,"/tmp_postQC.txt")) )
  file.remove(paste0(results_folder,"/tmp_postQC.txt"))
if( file.exists(paste0(results_folder,"/tmp_postQCAdj.txt")) )
  file.remove(paste0(results_folder,"/tmp_postQCAdj.txt"))

## Remove duplicates, Exclude CpGs and adjust data (BN and FDR)
for ( i in 1:length(files) )
{
  
  # Prepare output subfolder for cohort-model results (create if not exists)
  if(!dir.exists(file.path(getwd(), results_folder, prefixes[i] )))
    suppressWarnings(dir.create(file.path(getwd(), results_folder, prefixes[i])))
  
  # Creates an empty file to resume all data if an old file exist  removes
  # the file and creates a new one
  fResumeName <- paste0( file.path(getwd(), results_folder, prefixes[i]),"/",prefixes[i], "_descriptives.txt")
  if ( file.exists(fResumeName) ) {
    file.remove(fResumeName)
  }
  file.create(fResumeName)
  
  # Read data.
  cohort <- read.table(files[i], sep = "\t", header = TRUE, as.is = TRUE)
  
  # Ensure that probeID column exists properly, whether cpg is a column or a rowname
  if ("cpg" %in% colnames(cohort)) {
    colnames(cohort)[colnames(cohort) == "cpg"] <- "probeID"
  }
  else {
    cohort$probeID <- rownames(cohort)
  }
  
  if("beta" %in% colnames(cohort)) {
    colnames(cohort)[colnames(cohort) == "beta"] <- "BETA"
  }
  if("pval" %in% colnames(cohort)) {
    colnames(cohort)[colnames(cohort) == "pval"] <- "P_VAL"
  }
  
  print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))
  
  # Remove rows with NA from data
  cohort <- clean_NA_from_data(cohort)
  
  cohort <- descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], artype = artype[i], before = TRUE)
  
  
  test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
  
  # Remove cpGs with low representation
  # first, we test if colname_NforProbe and pcMissingSampes are defined
  if( !exists("colname_NforProbe") ) { colname_NforProbe <- NULL }
  if( !exists("pcMissingSamples") ) { pcMissingSamples <- NULL }
  
  cohort <- filterLowRepresentedCpGsinCohort(cohort, colname_NforProbe, pcMissingSamples, N[i], fileresume = fResumeName )
  
  # Exclude CpGs not meet conditions
  if("MASK_snp5_ethnic" %in% exclude ){
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = ethnic[i], filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
  } else {
    #..# if( !is.null(exclude) && exclude!='') {
    cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = "", filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
    #..# }
  }
  
  # Descriptives - After CpGs deletion #
  descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], artype = artype[i], before = FALSE )
  
  # Adjust data by Bonferroni and FDR
  cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, fResumeName, N[i]  )
  
  # Write QC complete data to external file
  write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))
  
  ## Visualization - Plots
  #. Problems in some workstations and servers.# rasterpdf::raster_pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'), res = 300)
  #..# Problems in some cases --> Get png plots : pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'))
  
  # Distribution plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_SE_plot.png'), type="cairo")
  plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
  dev.off()
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_pvals_plot.png'), type="cairo")
  plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')
  dev.off()
  
  # QQ plot
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_QQ_plot.png'), type="cairo")
  qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))
  dev.off()
  
  # Volcano plot.
  png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_Volcano_plot.png'), type="cairo")
  plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )
  dev.off()
  
  # Add mandatory data for precisionplot
  medianSE[i] <-  median(cohort$SE)
  value_N[i] <- N[i]
  cohort_label[i] <-  prefixes[i]
  
  # if n is defined for dichotomic condition :
  if(length(n) == length(N))  value_n[i] <- n[i]
  
  # Store data for Beta Box-Plot
  if( i == 1)
    betas.data <- list()
  betas.data[[prefixes[i]]] <- cohort[,"BETA"]
  
}
#


# Create QC Summary
if( file.exists(paste0(results_folder,"/tmp_pretQC.txt")) & file.exists(paste0(results_folder,"/tmp_postQC.txt")) & file.exists(paste0(results_folder,"/tmp_postQCAdj.txt") )) {
  preQC <- read.table (file = paste0(results_folder,"/tmp_pretQC.txt"), header = TRUE, sep = "\t")
  postQC <- read.table (file = paste0(results_folder,"/tmp_postQC.txt"), header = TRUE, sep = "\t")
  postQCAdj <- read.csv(file = paste0(results_folder,"/tmp_postQCAdj.txt"), header = TRUE, sep = "\t")
  
  write.table( cbind(prefixes,ethnic, postQC[,1:2], preQC, postQC[,3:length(postQC)], postQCAdj), file = paste0(results_folder,"/Summary_QCs_thesis.txt" ),
               row.names = FALSE, col.names = TRUE, sep = "\t")
  
  do.call(file.remove, list(list.files(results_folder, full.names = TRUE, pattern = "tmp_*")))
}


# Data for Precision Plot
precplot.data <- cbind.data.frame( SE = medianSE, invSE = (1/medianSE), N = value_N, sqrt_N = sqrt(N), cohort = cohort_label )
cols.numeric <- c("SE","invSE", "N", "sqrt_N")
precplot.data[cols.numeric] <- sapply(precplot.data[cols.numeric],as.numeric)

if(length(n) == length(N)){
  precplot.data.n <- cbind.data.frame( SE = medianSE, invSE = (1/medianSE), N = value_n, sqrt_N = sqrt(n), cohort = cohort_label )
  precplot.data.n[cols.numeric] <- sapply(precplot.data.n[cols.numeric],as.numeric)
}

# BoxPlot with Betas in all Models and cohorts
plot_betas_boxplot(betas.data, paste(results_folder, 'BETAS_BoxPlot.png', sep="/"))

#Boxplot per compound
patterns <- c("MEPA", "PRPA", "BPA", "TCS", "BP3", "ETPA")

#boxplots per model A B C
patterns <- c("Model1A", "Model1B", "Model1C")

#R
boxplot_folder <- "Z:/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis/Boxplots"
dir.create(boxplot_folder, showWarnings = FALSE)
setwd("Z:/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis/Boxplots")

#Mobaxterm
setwd("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis") 
boxplot_folder <- "/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis/Model_1/Boxplots_M1"
dir.create(boxplot_folder, showWarnings = FALSE)
setwd("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis/Model_1/Boxplots_M1")


for (pattern in patterns){
  grouped_data <- list()
  for (name in names(betas.data)){
    if (grepl(pattern, name)){
      grouped_data[[name]] <- betas.data[[name]]
      
    }
  }
  #moceh_names <- grep("MOCEH", names(grouped_data), value = TRUE)
  #grouped_data[moceh_names] <- NULL
  #option to leave one out
  #if (pattern == "BPA"){
  #moceh_names <- grep("MOCEH", names(grouped_data), value = TRUE)
  #grouped_data[moceh_names] <- NULL
  #}
  if (length(grouped_data) > 0){
    #open png
    file_path <- file.path(boxplot_folder, paste("Boxplot_", pattern, "_.png"))
    png(file_path, type = "cairo", width = 800, height = 600)
    par(mar = c(10, 5, 4, 2))  
    #find ylim 
    unlist_beta <- unlist(grouped_data)
    #if want scaled by outliers - remember to change in ylim= below
    #global_ylim <- range(unlist_beta, na.rm = TRUE)
    #if want scaled by IQR - remember to change in ylim= below 
    Q1 <- quantile(unlist_beta, 0.25, na.rm=TRUE)
    Q3 <- quantile(unlist_beta, 0.75, na.rm=TRUE)
    IQR_value <- Q3 - Q1 
    lower_bound <- Q1 - 5.5* IQR_value
    upper_bound <- Q3 + 5.5* IQR_value
    IQR_ylim <- c(lower_bound, upper_bound)
    #make boxplots
    boxplot(grouped_data, 
            main = paste("Boxplot for", pattern),  
            names = names(grouped_data),  # Label for the x-axis
            ylim = IQR_ylim,
            ylab = "Values",   # Label for the y-axis
            col = "lightblue", # Color for the boxplot
            outline = FALSE,   # Disable outliers (optional)
            range = 0, 
            xaxt= "n")
    axis(1, at = 1:length(names(grouped_data)), labels = names(grouped_data), las = 2, srt = 45)
    dev.off()
    cat("Saved combined boxplot for matched group", pattern, "\n")
  }
  
}

##  Post model analysis - precision plots  ## 
#R
precisionplots_folder <- "Z:/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis/Precision_Plots"
dir.create(precisionplots_folder, showWarnings = FALSE)
setwd("Z:/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis/Precision_Plots")

#Mobaxterm
setwd("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis") 
precisionplots_folder <- "/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis/Precision_Plots"
dir.create(precisionplots_folder, showWarnings = FALSE)
setwd("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis/Precision_Plots")


if ( length(files) > 1)
{
  # Precision_Plot(N) main 
  plot_precisionp(precplot.data, paste(precisionplots_folder,  "precision_SE_N.png",sep='/'), main = "Precision Plot - 1/median(SE) vs sqrt(n)")
  
  #divided by compound
  
  for (pattern in patterns){
    grouped_data <- data.frame()
    for (i in 1:nrow(precplot.data)){
      if (grepl(pattern, precplot.data$cohort[i])){
        grouped_data <- grouped_data <- rbind(grouped_data, data.frame(SE = precplot.data$SE[i],invSE = precplot.data$invSE[i],N = precplot.data$N[i], sqrt_N = precplot.data$sqrt_N[i],cohort = precplot.data$cohort[i]))
      }
    } 
    p <- ggplot2::ggplot( data = grouped_data, mapping = ggplot2::aes( x = round(as.numeric(sqrt_N),2),
                                                                       y = round(as.numeric(invSE),2) ) ) +
      ggplot2::theme_bw() +
      ggplot2::geom_point( size = 2, ggplot2::aes( colour = cohort ) ) +
      #geom_line( aes( group = cohort , colour = cohort ) ) +
      ggplot2::ggtitle( "Precision Plot - 1/median(SE) vs sqrt(n)" ) +
      ggplot2::theme( legend.position = "bottom",
                      legend.text = ggplot2::element_text(size=5),
                      legend.title = ggplot2::element_blank()) +
      ggplot2::labs( x = "sqrt(n)",
                     y = "inv SE") +
      ggplot2::scale_shape_manual( values = 0:7 ) +
      ggplot2::scale_x_continuous(expand = c(0.05, 0.05), limits = c(min(grouped_data$sqrt_N) - 0.1, max(grouped_data$sqrt_N) + 0.1)) +
      ggplot2::scale_y_continuous(expand = c(0.05, 0.05), limits = c(min(grouped_data$invSE) - 0.1, max(grouped_data$invSE) + 0.1))
    filename <- (paste(precisionplots_folder, paste(pattern, "_precision1_SE.png"), sep='/'))
    if(!is.null(filename)) {
      png(filename, type = "cairo")
      print(p)
      cat("Saved combined boxplot for matched group", pattern, "\n")
      dev.off()}
    
  }
  
  # Precision_Plot(n) main
  if(length(n) == length(N))
    plot_precisionp(precplot.data.n, paste(precisionplots_folder,  "precision1_SE_n.png", sep='/'), main = "Subgroup Precision Plot -  1/median(SE) vs sqrt(n)")
  
  # Venn_Diagrams()
  #for (i in 1:length(venn_diagrams))
  #plot_venndiagram(venn_diagrams[[i]], qcpath = results_folder, plotpath = results_folder, bn='padj.bonf', fdr='padj.fdr')
  
}

# Remove duplicates
# cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i], '/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )

-------------------------------------------------------------------------------------------------
  #useful table
  
setwd("Z:/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis")
setwd("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis/Model_1/Individual_results/") 

QC_sum <- read.table("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis/Model_1/Individual_results/Summary_QCs_thesis.txt", sep = "\t", header = TRUE)
View(QC_sum)
QC_sum_condensed <- QC_sum[ , c("prefixes", "NSamples", "Initial_CpGs", "N.Cpgs.after.QC", "Coef_min.","Coef_Mean", "Coef_max","SE_min.", "SE_Mean", "SE_max", "pval_min.", "pval_Mean", "lambda", "N_Sig_Nominal", "N_Sig_FDR", "N_Sig_BN")]
View(QC_sum_condensed)


QC_sum_sorted_by_compound <- QC_sum_condensed[order(grepl("MEPA$", QC_sum_condensed$prefixes), 
                                                    grepl("PRPA$", QC_sum_condensed$prefixes),
                                                    grepl("BPA", QC_sum_condensed$prefixes),
                                                    grepl("TCS$", QC_sum_condensed$prefixes),
                                                    grepl("BP3$", QC_sum_condensed$prefixes),
                                                    grepl("ETPA$", QC_sum_condensed$prefixes), 
                                                    QC_sum_condensed$prefixes),]
View(QC_sum_sorted_by_compound)


install.packages("writexl")
library(writexl)
write.csv(QC_sum_sorted_by_compound, file = "QC_sum_phenols_by_compound_final", row.names = FALSE )
write_xlsx(QC_sum_sorted_by_compound, "QC_phenol_excel_final.xlsx")



