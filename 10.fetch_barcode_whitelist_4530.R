library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4530"
sample="4530"
treatment="tumor"
figDir=paste0(homeDir,"project_results/figures/",sample,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))
rObjectsDir=paste0(homeDir,"project_results/Robjects/")
system(paste0("mkdir -p ",figDir))

dataH <- paste0("/share/ScratchGeneral/sunwu/Chromium_10X/count_data_2/CID4530/count_CID4520_GRCh38_mm10/outs/raw_gene_bc_matrices/GRCh38")
dataM <- paste0("/share/ScratchGeneral/sunwu/Chromium_10X/count_data_2/CID4530/count_CID4520_GRCh38_mm10/outs/raw_gene_bc_matrices/mm10")

seuratH <- Read10X(dataH)
seuratM <- Read10X(dataM)
seuratX<-rbind(seuratH,seuratM)
seuratMat<-as.matrix(seuratX[,colSums(seuratX)>500])
#this produces 6931 cells (cutoff 500)

write(colnames(seuratMat),file="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530/CITE_fastqs/CID4530-N_whitelist.txt")


dataH <- paste0("/share/ScratchGeneral/sunwu/Chromium_10X/count_data_2/CID4530-N/count_CID4520-N_GRCh38_mm10/outs/raw_gene_bc_matrices/GRCh38")
dataM <- paste0("/share/ScratchGeneral/sunwu/Chromium_10X/count_data_2/CID4530-N/count_CID4520-N_GRCh38_mm10/outs/raw_gene_bc_matrices/mm10")

seuratH <- Read10X(dataH)
seuratM <- Read10X(dataM)
seuratX<-rbind(seuratH,seuratM)
seuratMat<-as.matrix(seuratX[,colSums(seuratX)>500])
#this produces 6931 cells (cutoff 500)
write(colnames(seuratMat),file="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530/CITE_fastqs/CID4530_whitelist.txt")
