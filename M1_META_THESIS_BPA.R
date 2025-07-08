## ######################################### ##
##  Meta-Analysis to use with EASIER package ##
## ######################################### ##
#Model1 - BPA

# load package
library(EASIER)


########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to metaanalysis folder
setwd("Z:/analyses/ATH_EWAS_phtalates/")
setwd("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/")

# Files used in QC, needed in meta-analysis to plot ForestPlot
files_GenR_0y <- list.files(path = "./data/GENR/results/phenols/model1", 
                            pattern = "PACE_phthalates_GenR_EWAS_Model(A|B|C)_BPAugg_average_imputed_log2_.*", 
                            full.names = TRUE)
files_GenR_0y # 3 files


files_NELA <- list.files(path = "./data/NELA/results/phenols/model1", 
                         pattern = "PACE_phthalates_NELA_EWAS_Model(A|B|C)_(ad.*)_(BPA)(_media)?_S24_.*", 
                         full.names = TRUE)
files_NELA #3 files

files <- c(files_GenR_0y, files_NELA)
files #6 total files 

# Define a new directory for the modified files
output_dir <- "./data/renamed_files_forSFTHESIS_BPA"

# Create it if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through the files
renamed_files <- c()

for (file in files) {
  # Read original file
  cohort_x <- read.delim(file, header = TRUE)
  
  # Rename column if it exists
  if ("pval" %in% colnames(cohort_x)) {
    colnames(cohort_x)[colnames(cohort_x) == "pval"] <- "P_VAL"
  }
  
  # Define new file path in the output directory
  new_file <- file.path(output_dir, basename(file))
  
  # Save modified file
  write.table(cohort_x, file = new_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Keep track of new file paths
  renamed_files <- c(renamed_files, new_file)
}



# Prefixes for each file
prefixes <- c(
  'Model1A_GenR_BPA','Model1B_GenR_BPA', 'Model1C_GenR_BPA',
   'Model1A_NELA_BPA',
   'Model1B_NELA_BPA',
   'Model1C_NELA_BPA'
)
length(prefixes)


# Samples in original files used in QC
N <- c(
  276, 276, 276, #GenR0y, changed already
  314, 314, 314 #NELA
)
length(N)

# Define data for each meta-analysis
metafiles <- list(
  'Meta1A' = c('Model1A_GenR_BPA',  
              'Model1A_NELA_BPA' ),
  
  'Meta1B' = c('Model1B_GenR_BPA',  
              'Model1B_NELA_BPA'),
  
  'Meta1C' = c( 'Model1C_GenR_BPA',  
                'Model1C_NELA_BPA'))


# Array type, used in each meta-analysis : EPIC, 450K or MIX (to be used when in meta-analyses we have 450K and EPIC arrays)
artype <- c("MIX", "MIX", "MIX")


# Define maximum percent missing for each CpG
#     if pcenMissin = 0 only runs meta-analysis with all data
pcentMissing <- 0.8 # CpGs with precense lower than pcentMissing after GWAS meta-analysis will be deleted from the study.


# Paths with QCResults and path to store GWAMA results

# Paths with QCResults and path to store GWAMA results
results_folder <- 'metaanalysis/QC_Results_Phenols_Sfraley_thesis/Model_1/Individual_results'
results_gwama <- 'metaanalysis/GWAMA_Results_Phenols_Sfraley_thesis/Model_1_BPA'


# Venn diagrams ==> IMPORTANT : maximum 5 meta-analysis in each venn diagram
venndiag_threshold <- 0.05
venn_diagrams <- list(
  c("Meta1A", "Meta1B", "Meta1C" ),
  c("Meta1A_Filtr", "Meta1B_Filtr", "Meta1C_Filtr" )
)

########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########

# GWAMA binary path  (GWAMA IsGlobal Server installation)
# gwama.dir <- paste0(Sys.getenv("HOME"), "/data/software/GWAMA/")
gwama.dir <- "/PROJECTES/PUBLICDATA/software/EASIER/GWAMA/"


## Create directory for GWAMA configuration files and GWAMA_Results
if(!dir.exists(file.path(paste(results_gwama, "GWAMA", sep="/") )))
  suppressWarnings(dir.create(file.path( paste(results_gwama, "GWAMA", sep="/")), recursive = TRUE))

## Create directory for GWAMA_Results
outputfolder <- paste0(results_gwama, "/GWAMA_Results")
if(!dir.exists(file.path( outputfolder )))
  suppressWarnings(dir.create(file.path( outputfolder), recursive = TRUE))


# Create hapmap files for the different artypes that we cab use (450K and EPIC)
# map file is used in Manhattan plots
hapmapfile_450K <- paste(results_gwama,"GWAMA", "hapmap_450K.map" ,sep = "/")
generate_hapmap_file("450K", hapmapfile_450K)
hapmapfile_EPIC <- paste(results_gwama,"GWAMA", "hapmap_EPIC.map" ,sep = "/")
generate_hapmap_file("EPIC", hapmapfile_EPIC)
hapmapfile_MIX <- paste(results_gwama,"GWAMA", "hapmap_MIX.map" ,sep = "/")
generate_hapmap_file("MIX", hapmapfile_MIX)


for( metf in 1:length(metafiles))
{
  
  list.lowCpGs <- NULL
  
  # Create folder for a meta-analysis in GWAMA folder, here we store the GWAMA input files for each meta-analysis,
  # We create one for complete meta-analysis
  if(!dir.exists(file.path( paste(results_gwama,"GWAMA", names(metafiles)[metf] ,sep="/") )))
    suppressWarnings(dir.create(file.path( paste(results_gwama,"GWAMA", names(metafiles)[metf], sep="/"))))
  # We create another for meta-analysis without filtered CpGs with low percentage (sufix _Filtr)
  if(!dir.exists(file.path( paste0(results_gwama,"/GWAMA/", names(metafiles)[metf],"_Filtr") )))
    suppressWarnings(dir.create(file.path( paste0(results_gwama,"/GWAMA/", names(metafiles)[metf],"_Filtr"))))
  
  # GWAMA File name base
  inputfolder <- paste0(results_gwama,"/GWAMA/",  names(metafiles)[metf])
  
  modelfiles <- unlist(metafiles[metf])
  
  runs <- c('Normal', 'lowcpgs') # Execution with all CpGs and without filtered CpGs
  lowCpGs = FALSE;
  outputfiles <- list()
  
  outputgwama <- paste(outputfolder,names(metafiles)[metf],sep = '/')
  
  for(j in 1:length(runs))
  {
    if(runs[j]=='lowcpgs') {
      lowCpGs = TRUE
      # Get low presence CpGs in order to exclude this from the new meta-analysis
      list.lowCpGs <- get_low_presence_CpGs(outputfiles[[j-1]], pcentMissing)
      inputfolder <- paste0(results_gwama,"/GWAMA/",  names(metafiles)[metf], "_Filtr")
      outputgwama <- paste0(outputgwama,"_Filtr")
    }
    
    # Create GWAMA files for each file in meta-analysis and execute GWAMA
    for ( i in 1:length(modelfiles) )
      create_GWAMA_files(file.path(results_folder,modelfiles[i]),  modelfiles[i], inputfolder, N[which(prefixes==modelfiles[i])], list.lowCpGs )
    
    
    # Get hapmapfile attending to current metaanalysis artype
    hapmapfile <- hapmapfile_450K
    if(artype[metf]=='EPIC'){
      hapmapfile <- hapmapfile_EPIC
    } else if(artype[metf]=='MIX'){
      hapmapfile <- hapmapfile_MIX
    }
    
    #.Original.#outputfiles[[runs[j]]] <- execute_GWAMA_MetaAnalysis(prefixgwama, names(metafiles)[metf])
    outputfiles[[runs[j]]] <- run_GWAMA_MetaAnalysis(inputfolder, outputgwama, names(metafiles)[metf], gwama.dir, hapmapfile)
    
    # Post-metha-analysis QC --- >>> adds BN and FDR adjustment
    dataPost <- get_descriptives_postGWAMA(outputgwama, outputfiles[[runs[j]]], modelfiles, names(metafiles)[metf], artype[metf], N[which(prefixes %in% modelfiles)] )
    
    # Forest-Plot
    plot_ForestPlot( dataPost, metafiles[[metf]], runs[j], inputfolder, names(metafiles)[metf], renamed_files, outputgwama, nsignificatives = 30  )
    
  }
  
}


# Venn_Diagrams for for meta-analysis with fixed effects

for (i in 1:length(venn_diagrams)){
  if(length(venn_diagrams[[i]])>1){
    plot_venndiagram(venn_diagrams[[i]], qcpath = outputfolder, plotpath =  paste0(results_gwama, "/GWAMA_Results"), pattern = '_Fixed_Modif.out',bn='Bonferroni', fdr='FDR', venndiag_threshold)
  }
}

if(dir.exists(file.path( paste(results_gwama, "GWAMA", sep="/") )))
  unlink(file.path(results_gwama, "GWAMA"), recursive=TRUE)

