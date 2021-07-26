#### Installation ####

#Step 1: Install package dependencies 
metanr_packages <- function(){
  
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn","httr","qs")
  
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs)!=0){
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

#Step 2: Install the package

# Step 1: Install devtools
install.packages("devtools")
library(devtools)

# Step 2: Install MetaboAnalystR with documentation
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T)

devtools::install_github("xia-lab/OptiLCMS", build = TRUE, build_vignettes = FALSE, build_manual =T)
library(OptiLCMS)

# Load the MetaboAnalystR package
library("MetaboAnalystR")

## Setting the data depositing folder
data_folder_Sample <- "~/Data_IBD"
data_folder_QC <- "~/QC_IBD"  

# Use Google API for data downloading here. 
# Please "install.packages('googledrive')" and "install.packages('httpuv')"first.
library(googledrive);
temp <- tempfile(fileext = ".zip")
# Please authorize your google account to access the data
dl <- drive_download(
  as_id("10DBpPEWy2cZyXvmlLOIfpwqQYwFplKYK"), path = temp, overwrite = TRUE)
# Setting your own date file folder
out <- unzip(temp, exdir = data_folder_QC)
# Date files for parameters optimization are deposited below
out

# Inspect the MS data via a 3D image. "res" are used to specify the resolution for the MS data.
PerformDataInspect(data_folder_QC,res = 50)

# QC samples are trimmed with "ssm_trim" strategy. 
# The percentage of RT dimension retained is set as 20%.
raw_data <- PerformDataTrimming(data_folder_QC,rt.idx = 0.2)

# Initial platform specific parameters
param_initial <- SetPeakParam(platform = "UPLC-Q/E") 

param_optimized <- PerformParamsOptimization(raw_data, param = param_initial, ncore = 8) 

# Import raw MS data. The "SetPlotParam" parameters can be used to determine plot or not.
rawData <- ImportRawMSData(data_folder_Sample,plotSettings = SetPlotParam(Plot=F))
