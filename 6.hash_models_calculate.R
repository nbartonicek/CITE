library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4497"
sample="CID4497-PDX_Hash"
treatment="original"
treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sample,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))
rObjectsDir=paste0(homeDir,"project_results/Robjects/")
system(paste0("mkdir -p ",figDir))

load("../project_results/Robjects/Human_HASHSEQ_S1_hd2_F_clean.Rdata")
load(="../project_results/Robjects/Human_HASHSEQ_S1_hd2_F_seuratwithMito.Rdata")
