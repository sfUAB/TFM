## ############ ##
##  Enrichment  ##
## ############ ##
#Sylvie Fraley - Model 1 THESIS BPA

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


# Load package
require(EASIER)


## Develop test working directory
# setwd("~/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/EASIER")
# setwd("/Users/mailos/tmp/proves")

########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to enrichment folder
setwd("Z:\analyses\ATH_EWAS_phtalates\metaanalysis\Enrichment_Sfraley_thesis")
setwd("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/Enrichment_Sfraley_thesis")


# Files with CpG data to enrich may be a CpGs list or annotated GWAMA output
FilesToEnrich <- c("/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/GWAMA_Results_Phenols_Sfraley_thesis/Model_1_BPA/GWAMA_Results/Meta1C_Filtr/Meta1C_Fixed_Modif.out",
                   "/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/GWAMA_Results_Phenols_Sfraley_thesis/Model_1_BP3/GWAMA_Results/Meta1C_Filtr/Meta1C_Fixed_Modif.out",
                    "/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/GWAMA_Results_Phenols_Sfraley_thesis/Model_1_MEPA/GWAMA_Results/Meta1C_Filtr/Meta1C_Fixed_Modif.out",
                   "/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/GWAMA_Results_Phenols_Sfraley_thesis/Model_1_PRPA/GWAMA_Results/Meta1C_Filtr/Meta1C_Fixed_Modif.out"
)


# Values for adjustment

BN <-  TRUE    # change to alternate one? 
FDR <- 0.7     # significance level for adjustment, if NA FDR is not used
pvalue <- 0.00001 # significance level for p-value, if NA p-value is not used, suggestive ones

# Array type, used : EPIC or 450K
# this data is defined for each file to analyse
artype <- c('450K', '450K', '450K') 

# Result paths definition for QC, Meta-Analysis and Enrichment
results_folder <- "/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/QC_Results_Phenols_Sfraley_thesis/Model_1/Individual_results"
results_gwama <- "/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_phtalates/metaanalysis/GWAMA_Results_Phenols_Sfraley_thesis/" #not sure
results_enrich <- "." #correct if at current wd

# Enrichment type :  'BLOOD' or 'PLACENTA'
#     if enrichtype <- 'BLOOD' => enrichment with :
#                          Cromatine States : BLOOD (crom15)
#                          (To be implemented in future) Partially Methylated Domains (PMD) for Blood
#     if enrichtype <- 'PLACENTA' => enrichment with:
#                          Cromatine States : PLACENTA (FP_15) optionally (FP_18)
#                          Partially Methylated Domains (PMD) for Placenta
#     if enrichtype is different from 'BLOOD' and 'PLACENTA' we only get the missMethyl and MSigDB enrichment and the Unique genes list.
enrichtype <- 'BLOOD'

# Cromatine States Placenta Enrichment FP_18
# if enrichFP18 = TRUE the enrichment is performed wit FP_15 and FP_18
enrichFP18 <- FALSE

# Test to be used : 'Fisher' or 'Hypergeometric' if testdata is different no test will be performed
testdata <- 'Fisher'

# Perform eQTM enrichment
bEQTM <- TRUE

########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########


## Check if we have any files to enrich and if these files exists
if (length(FilesToEnrich)>=1 & FilesToEnrich[1]!='') {
  for ( i in 1:length(FilesToEnrich))
    if (!file.exists(FilesToEnrich[i])) stop(paste0('File ',FilesToEnrich[i],' does not exsits, please check file' ))
}

## Check variables

if( ! toupper(enrichtype) %in% c('PLACENTA','BLOOD') )
  warning('Only enrichment with MyssMethyl and MSigDB will be done')

if( ! tolower(testdata) %in% c('fisher','hypergeometric') )
  warning('Wrong value for testdata variable, values must be "Fisher" or "Hypergeometric". No test will be performed ')



# Convert relative paths to absolute paths for FilesToEnrich
FilesToEnrich <- unlist(sapply(FilesToEnrich, function(file) { if(substr(file,1,1)!='.' & substr(file,1,1)!='/' & substr(file, 2, 2) != ':') file <- paste0('./',file) else file }))
FilesToEnrich <- sapply(FilesToEnrich, tools::file_path_as_absolute)

if(results_enrich!='.'){
  outputfolder <- file.path(getwd(), results_enrich )
}else{
  outputfolder <- file.path(getwd() )}


# Create dir to put results from enrichment
if(!dir.exists(outputfolder))
  suppressWarnings(dir.create(outputfolder))

setwd( outputfolder)

# Get which data we have to enrich
if (length(FilesToEnrich)>=1 & FilesToEnrich[1]!='')
{
  
  for (i in 1:length(FilesToEnrich)) {
    
    # Enrich all CpGs
    allCpGs <- FALSE
    
    # Get data
    data <- NULL
    data <- read.table(FilesToEnrich[i], header = TRUE, sep = "", dec = ".", stringsAsFactors = FALSE)
    
    # Is a CpG list only ? then read without headers and annotate data
    if(nrow(data) <= 1 || ncol(data) <= 1 || !any(c("FDR", "BN","p.value", "Bonferroni") %in% colnames(data)) ) {
      data <- read.table(FilesToEnrich[i], dec = ".") # Avoid header
      data <- as.vector(t(data))
      data <- get_annotattions(data, artype[i], FilesToEnrich[i], outputfolder )
      allCpGs <- TRUE
      data$chromosome <- substr(data$chr,4,length(data$chr))
      data$rs_number <- data$CpGs
    }else {
      if(! "rs_number" %in% colnames(data)) {
        if("CpGs" %in% colnames(data)) {
          data$rs_number = data$CpGs
        }else if("CpGId" %in% colnames(data)) {
          data$rs_number = data$CpGId
        }else {
          stop("Data must contain rs_number, CpGs or CpGId column with CpGs Ids")
        }
      }
    }
    
    ## -- Functional Enrichmnet
    ## ------------------------
    
    # Enrichment with missMethyl - GO and KEGG --> Writes results to outputfolder
    miss_enrich <- missMethyl_enrichment(data, outputfolder, FilesToEnrich[i], artype[i], BN, FDR, pvalue, allCpGs, plots = FALSE )
    
    # get unique genes from data
    geneUniv <- lapply( lapply(miss_enrich[grepl("signif", names(miss_enrich))], function(cpgs) { data[which(as.character(data$CpGs) %in% cpgs),]$UCSC_RefGene_Name}), getUniqueGenes)
    
    
    ## -- Online Tools
    
    # Enrichment with ConsensusPathDB
    #     - Consensus path http://cpdb.molgen.mpg.de/ (gene-set analysis – over-representation analysis)
    
    # Available FSet types :
    # 1 P     manually curated pathways from pathway databases
    # 2 N     interaction network neighborhood-based functional sets
    # 3 G2    Gene Ontology-based sets, GO level 2
    # 4 G3    Gene Ontology-based sets, GO level 3
    # 5 G4    Gene Ontology-based sets, GO level 4
    # 6 G5    Gene Ontology-based sets, GO level 5
    # 7 C     protein complex-based sets
    
    acFSet <- c('C', 'P', 'G2', 'G3')
    acType <- 'entrez-gene'
    
    # Get Enrichment
    CPDB_enrich <- lapply(names(geneUniv), function( data, accFSet, genes ) {
      print(data)
      lapply(accFSet,
             get_consensusPdb_OverRepresentation,
             entityType='genes',
             accNumbers=na.omit(as.character(eval(parse(text = paste0("genes$",data))))),
             accType=acType,
             outputdir = "ConsensusPathDB",
             outputfile = gsub(".", "_", data, fixed=TRUE) )},
      accFSet = acFSet, genes = geneUniv)
    
    names(CPDB_enrich) <- names(geneUniv)
    
    
    ## -- Molecular Enrichmnet
    ## -----------------------
    
    # Molecular Signatures Database enrichment
    msd_enrich <- MSigDB_enrichment(data, outputfolder, FilesToEnrich[i], artype[i], BN, FDR, pvalue, allCpGs)
    
    
    if("FDR" %in% colnames(data) | "Bonferroni" %in% colnames(data) |  "p.val" %in% colnames(data) | "p.value" %in% colnames(data))
    {
      
      ## -- Prepare data
      ## ---------------
      
      # Classify by Hyper and Hypo methylated
      data$meth_state <- getHyperHypo(data$beta) # Classify methylation into Hyper and Hypo
      
      
      if("FDR" %in% colnames(data) & !is.na(FDR) )
      {
        # Add column bFDR to data for that CpGs that accomplish with FDR
        data$bFDR <- getBinaryClassificationYesNo(data$FDR, "<", FDR) # Classify fdr into "yes" and no taking into account FDR significance level
        
        # CpGs FDR and Hyper and Hypo respectively
        FDR_Hyper <- ifelse(data$bFDR == 'yes' & data$meth_state=='Hyper', "yes", "no")
        FDR_Hypo <- ifelse(data$bFDR == 'yes' & data$meth_state=='Hypo', "yes", "no")
        FDR_Hyper_Hypo <- ifelse(data$bFDR == 'yes' & data$meth_state=='Hypo', "Hypo-yes",
                                 ifelse(data$bFDR == 'no' & data$meth_state=='Hypo', "Hypo-no",
                                        ifelse(data$bFDR == 'yes' & data$meth_state=='Hyper', "Hyper-yes", "Hyper-no" ) ) )
      }
      
      if("Bonferroni" %in% colnames(data) & BN==TRUE)
      {
        # CpGs Bonferroni and Hyper and Hypo respectively
        BN_Hyper <- ifelse(data$Bonferroni == 'yes' & data$meth_state=='Hyper', "yes", "no")
        BN_Hypo <- ifelse(data$Bonferroni == 'yes' & data$meth_state=='Hypo', "yes", "no")
        BN_Hyper_Hypo <- ifelse(data$Bonferroni == 'yes' & data$meth_state=='Hypo', "Hypo-yes",
                                ifelse(data$Bonferroni == 'no' & data$meth_state=='Hypo', "Hypo-no",
                                       ifelse(data$Bonferroni == 'yes' & data$meth_state=='Hyper', "Hyper-yes", "Hyper-no" ) ) )
      }
      
      
      if(("p.val" %in% colnames(data) | "p.value" %in% colnames(data)) & !is.na(pvalue) )
      {
        
        # Add column bpval to data for that CpGs that accomplish with FDR
        data$bpval <- getBinaryClassificationYesNo(data$p.value, "<", pvalue) # Classify fdr into "yes" and no taking into account FDR significance level
        
        # CpGs FDR and Hyper and Hypo respectively
        pval_Hyper <- ifelse(data$bpval == 'yes' & data$meth_state=='Hyper', "yes", "no")
        pval_Hypo <- ifelse(data$bpval == 'yes' & data$meth_state=='Hypo', "yes", "no")
        pval_Hyper_Hypo <- ifelse(data$bpval == 'yes' & data$meth_state=='Hypo', "Hypo-yes",
                                  ifelse(data$bpval == 'no' & data$meth_state=='Hypo', "Hypo-no",
                                         ifelse(data$bpval == 'yes' & data$meth_state=='Hyper', "Hyper-yes", "Hyper-no" ) ) )
      }
      
      
      ## TODO: Simplify this code with only one function x option (fisher - Geometric // BN - FDR )
      
      # For FDR
      if("FDR" %in% colnames(data) & !is.na(FDR) )
      {
        
        ## --  CpG Gene position
        ## ---------------------
        
        # Get descriptives
        get_descriptives_GenePosition(data$UCSC_RefGene_Group, data$bFDR , "FDR", outputdir = "GenePosition/Fisher_FDR_Desc", outputfile = FilesToEnrich[i])
        
        if( tolower(testdata) =='fisher') {
          ## --  Fisher Test - Gene position - FDR, FDR_hyper and FDR_hypo
          GenePosition <- getAllFisherTest(data$bFDR, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_FDR", outputfile = FilesToEnrich[i], plots = TRUE )
          GenePosition_hyper <- getAllFisherTest(FDR_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_FDRHyper", outputfile = FilesToEnrich[i], plots = TRUE )
          GenePosition_hypo <- getAllFisherTest(FDR_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_FDRHypo", outputfile = FilesToEnrich[i], plots = TRUE )
          GenePosition_HyperHypo <- getAllFisherTest(FDR_Hyper_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_FDRHyperHypo", outputfile = FilesToEnrich[i], plots = TRUE )
        }
        else if ( tolower(testdata) =='hypergeometric') {
          ## --  HyperGeometric Test - Island relative position - FDR, FDR_hyper and FDR_hypo (for Depletion and Enrichment)
          GenePosition <- getAllHypergeometricTest(data$bFDR, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_FDR", outputfile = FilesToEnrich[i])
          GenePosition_hyper <- getAllHypergeometricTest(FDR_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_FDRHyper", outputfile = FilesToEnrich[i])
          GenePosition_hypo <- getAllHypergeometricTest(FDR_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_FDRHypo", outputfile = FilesToEnrich[i])
        }
        
        plot_TestResults_Collapsed(list(fdr = GenePosition, fdr_hypo = GenePosition_hypo, fdr_hyper = GenePosition_hyper),
                                   outputdir = "GenePosition", outputfile = FilesToEnrich[i], main = )
        
        ## --  CpG Island relative position
        ## --------------------------------
        
        # Get descriptives
        get_descriptives_RelativetoIsland(data$Relation_to_Island, data$bFDR , "FDR", outputdir = "RelativeToIsland/Fisher_FDR_RelativeToIsland", outputfile = FilesToEnrich[i])
        
        if( tolower(testdata) =='fisher') {
          ## --  Fisher Test - Position Relative to Island - FDR, FDR_hyper and FDR_hypo
          relative_island <- getAllFisherTest(data$bFDR, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_FDR", outputfile = FilesToEnrich[i], plots = TRUE )
          relative_island_hyper <- getAllFisherTest(FDR_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_FDRHyper", outputfile = FilesToEnrich[i], plots = TRUE )
          relative_island_hypo <- getAllFisherTest(FDR_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_FDRHypo", outputfile = FilesToEnrich[i], plots = TRUE )
          relative_island_hyperhypo <- getAllFisherTest(FDR_Hyper_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_FDRHyperHypo", outputfile = FilesToEnrich[i], plots = TRUE )
        }
        else if ( tolower(testdata) =='hypergeometric') {
          ## --  HyperGeometric Test - Gene position - FDR, FDR_hyper and FDR_hypo (for Depletion and Enrichment)
          relative_island <- getAllHypergeometricTest(data$bFDR, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_FDR", outputfile = FilesToEnrich[i])
          relative_island_hyper <- getAllHypergeometricTest(FDR_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_FDRHyper", outputfile = FilesToEnrich[i])
          relative_island_hypo <- getAllHypergeometricTest(FDR_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_FDRHypo", outputfile = FilesToEnrich[i])
        }
        plot_TestResults_Collapsed(list(bn = relative_island, bn_hypo = relative_island_hypo, bn_hyper = relative_island_hyper),
                                   outputdir = "RelativeToIsland", outputfile = FilesToEnrich[i], main = )
        
        
      }
      
      # For Bonferroni
      if("Bonferroni" %in% colnames(data) & BN==TRUE)
      {
        
        ## --  CpG Gene position
        ## ---------------------
        
        get_descriptives_GenePosition(data$UCSC_RefGene_Group, data$Bonferroni, "Bonferroni", outputdir = "GenePosition/Fisher_BN_Desc", outputfile = FilesToEnrich[i])
        
        # For BN
        if( tolower(testdata) =='fisher') {
          ## --  Fisher Test - Gene position - FDR, FDR_hyper and FDR_hypo
          GenePosition <- getAllFisherTest(data$Bonferroni, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_BN", outputfile = FilesToEnrich[i], plots = TRUE )
          GenePosition_hyper <- getAllFisherTest(BN_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_BNHyper", outputfile = FilesToEnrich[i], plots = TRUE )
          GenePosition_hypo <- getAllFisherTest(BN_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_BNHypo", outputfile = FilesToEnrich[i], plots = TRUE )
          GenePosition_hyperhypo <- getAllFisherTest(BN_Hyper_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_BNHyperHypo", outputfile = FilesToEnrich[i], plots = TRUE )
        }
        else if ( tolower(testdata) =='hypergeometric') {
          ## --  HyperGeometric Test - Island relative position - FDR, FDR_hyper and FDR_hypo (for Depletion and Enrichment)
          GenePosition <- getAllHypergeometricTest(data$Bonferroni, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_BN", outputfile = FilesToEnrich[i])
          GenePosition_hyper <- getAllHypergeometricTest(BN_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_BNHyper", outputfile = FilesToEnrich[i])
          GenePosition_hypo <- getAllHypergeometricTest(BN_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_BNHypo", outputfile = FilesToEnrich[i])
        }
        
        plot_TestResults_Collapsed(list(fdr = GenePosition, fdr_hypo = GenePosition_hypo, fdr_hyper = GenePosition_hyper),
                                   outputdir = "GenePosition", outputfile = FilesToEnrich[i], main = )
        
        ## --  CpG Island relative position
        ## --------------------------------
        
        # Get descriptives
        get_descriptives_RelativetoIsland(data$Relation_to_Island, data$Bonferroni, "Bonferroni", outputdir = "RelativeToIsland/Fisher_BN_RelativeToIsland", outputfile = FilesToEnrich[i])
        
        if( tolower(testdata) =='fisher') {
          ## --  Fisher Test - Position Relative to Island - FDR, FDR_hyper and FDR_hypo
          relative_island <- getAllFisherTest(data$Bonferroni, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_BN", outputfile = FilesToEnrich[i], plots = TRUE )
          relative_island_hyper <- getAllFisherTest(BN_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_BNHyper", outputfile = FilesToEnrich[i], plots = TRUE )
          relative_island_hypo <- getAllFisherTest(BN_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_BNHypo", outputfile = FilesToEnrich[i], plots = TRUE )
          relative_island_hyperhypo <- getAllFisherTest(BN_Hyper_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_BNHyperHypo", outputfile = FilesToEnrich[i], plots = TRUE )
        }
        else if ( tolower(testdata) =='hypergeometric') {
          ## --  HyperGeometric Test - Gene position - FDR, FDR_hyper and FDR_hypo (for Depletion and Enrichment)
          relative_island <- getAllHypergeometricTest(data$Bonferroni, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_BN", outputfile = FilesToEnrich[i])
          relative_island_hyper <- getAllHypergeometricTest(BN_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_BNHyper", outputfile = FilesToEnrich[i])
          relative_island_hypo <- getAllHypergeometricTest(BN_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_BNHypo", outputfile = FilesToEnrich[i])
        }
        plot_TestResults_Collapsed(list(bn = relative_island, bn_hypo = relative_island_hypo, bn_hyper = relative_island_hyper),
                                   outputdir = "RelativeToIsland", outputfile = FilesToEnrich[i], main = )
        
      }
      
      # For pvalue
      if( ("p.val" %in% colnames(data) | "p.value" %in% colnames(data)) & !is.na(pvalue) )
      {
        
        ## --  CpG Gene position
        ## ---------------------
        
        # Get descriptives
        get_descriptives_GenePosition(data$UCSC_RefGene_Group, data$bpval , "p.value", outputdir = "GenePosition/Fisher_pval_Desc", outputfile = FilesToEnrich[i])
        
        if( tolower(testdata) =='fisher') {
          ## --  Fisher Test - Gene position - pval, pval_hyper and pval_hypo
          GenePosition <- getAllFisherTest(data$bpval, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_pval", outputfile = FilesToEnrich[i], plots = TRUE )
          GenePosition_hyper <- getAllFisherTest(pval_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_pvalHyper", outputfile = FilesToEnrich[i], plots = TRUE )
          GenePosition_hypo <- getAllFisherTest(pval_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_pvalHypo", outputfile = FilesToEnrich[i], plots = TRUE )
          GenePosition_hyperhypo <- getAllFisherTest(pval_Hyper_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_pvalHyperHypo", outputfile = FilesToEnrich[i], plots = TRUE )
        }
        else if ( tolower(testdata) =='hypergeometric') {
          ## --  HyperGeometric Test - Island relative position - pval, pval_hyper and pval_hypo (for Depletion and Enrichment)
          GenePosition <- getAllHypergeometricTest(data$bpval, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_pval", outputfile = FilesToEnrich[i])
          GenePosition_hyper <- getAllHypergeometricTest(pval_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_pvalHyper", outputfile = FilesToEnrich[i])
          GenePosition_hypo <- getAllHypergeometricTest(pval_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_pvalHypo", outputfile = FilesToEnrich[i])
        }
        
        plot_TestResults_Collapsed(list(pval = GenePosition, pval_hypo = GenePosition_hypo, pval_hyper = GenePosition_hyper),
                                   outputdir = "GenePosition", outputfile = FilesToEnrich[i], main = )
        
        ## --  CpG Island relative position
        ## --------------------------------
        
        # Get descriptives
        get_descriptives_RelativetoIsland(data$Relation_to_Island, data$bpval , "p.value", outputdir = "RelativeToIsland/Fisher_pval_RelativeToIsland", outputfile = FilesToEnrich[i])
        
        if( tolower(testdata) =='fisher') {
          ## --  Fisher Test - Position Relative to Island - pval, pval_hyper and pval_hypo
          relative_island <- getAllFisherTest(data$bpval, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_pval", outputfile = FilesToEnrich[i], plots = TRUE )
          relative_island_hyper <- getAllFisherTest(pval_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_pvalHyper", outputfile = FilesToEnrich[i], plots = TRUE )
          relative_island_hypo <- getAllFisherTest(pval_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_pvalHypo", outputfile = FilesToEnrich[i], plots = TRUE )
          relative_island_hyperhypo <- getAllFisherTest(pval_Hyper_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_pvalHyperHypo", outputfile = FilesToEnrich[i], plots = TRUE )
        }
        else if ( tolower(testdata) =='hypergeometric') {
          ## --  HyperGeometric Test - Gene position - pval, pval_hyper and pval_hypo (for Depletion and Enrichment)
          relative_island <- getAllHypergeometricTest(data$bpval, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_pval", outputfile = FilesToEnrich[i])
          relative_island_hyper <- getAllHypergeometricTest(pval_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_pvalHyper", outputfile = FilesToEnrich[i])
          relative_island_hypo <- getAllHypergeometricTest(pval_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_pvalHypo", outputfile = FilesToEnrich[i])
        }
        plot_TestResults_Collapsed(list(bn = relative_island, bn_hypo = relative_island_hypo, bn_hyper = relative_island_hyper),
                                   outputdir = "RelativeToIsland", outputfile = FilesToEnrich[i], main = )
        
        
      }
      
      
    } else {
      
      ## -- Prepare data
      ## ---------------
      
      # * Create dataframe with non significative CpGs attending to artype ('EPIC' or '450K')
      #     - Get gene position tests
      #     - Get CpG Island relative position ests
      
      # Get
      unsignif_df <- get_annotation_unlisted_CpGs(data$rs_number, artype[i])
      
      data$signif <- 'yes'
      unsignif_df$signif <- 'no'
      
      data <- rbind(data[,which(colnames(data) %in% colnames(as.data.frame(unsignif_df)))], as.data.frame(unsignif_df) )
      data$rs_number <- data$Name
      
      
      ## TODO: Simplify this code with only one function x option (fisher - Geometric // BN - FDR )
      
      ## --  CpG Gene position
      ## ---------------------
      
      # Get descriptives
      get_descriptives_GenePosition(data$UCSC_RefGene_Group, data$signif , "CpGlist", outputdir = "GenePosition/Fisher_CpGlist_Desc", outputfile = FilesToEnrich[i])
      
      if( tolower(testdata) =='fisher') {
        GenePosition <- getAllFisherTest(data$signif, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_CpGlist", outputfile = FilesToEnrich[i], plots = TRUE )
      }else if ( tolower(testdata) =='hypergeometric') {
        GenePosition <- getAllHypergeometricTest(data$signif, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_CpGlist", outputfile = FilesToEnrich[i])
      }
      
      #..# plot_RelativetoIsland(GenePosition, outputdir = "GenePosition", outputfile = paste0("CpGlist_",FilesToEnrich[i]), main = )
      plot_GenePosition(GenePosition, outputdir = "GenePosition", outputfile = paste0("CpGlist_",FilesToEnrich[i]), main = )
      
      ## --  CpG Island relative position
      ## --------------------------------
      
      # Get descriptives
      get_descriptives_RelativetoIsland(data$Relation_to_Island, data$signif , "CpGlist", outputdir = "RelativeToIsland/Fisher_CpGlist_RelativeToIsland", outputfile = FilesToEnrich[i])
      
      if( tolower(testdata) =='fisher') {
        relative_island <- getAllFisherTest(data$signif, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_CpGlist", outputfile = FilesToEnrich[i], plots = TRUE )
      } else {
        relative_island <- getAllHypergeometricTest(data$signif, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_CpGlist", outputfile = FilesToEnrich[i])
      }
      
      plot_OR(relative_island, outputdir = "RelativeToIsland", outputfile = paste0("CpGlist_",FilesToEnrich[i]), main = )
      
      #..# plot_TestResults_Collapsed(list(relat = relative_island),
      #..#                            outputdir = "RelativeToIsland", outputfile = paste0("CpGlist_",FilesToEnrich[i]), main = )
      
    }
    
  
    # WRITE FINAL ENRICHMENT DATA
    write.table( crom_data, paste0( getwd(), "/",tools::file_path_sans_ext(basename(FilesToEnrich[i])),"_Enriched.csv" ) , quote=F, row.names=F, sep="\t")
  }
  
} else{
  print ("Error no data to enrich.")
}
