library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(DropletUtils)


homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="3921_3963_4066"
sample="3921_3963_4066"
treatment="original"
figDir=paste0(homeDir,"project_results/figures/",sample,"/")
system(paste0("mkdir -p ",figDir))
rObjectsDir=paste0(homeDir,"project_results/Robjects/")
system(paste0("mkdir -p ",figDir))

#emptyDrops


dataH <- paste0(homeDir,"raw_files/",sample,"/outs/raw_gene_bc_matrices/GRCh38")

seuratX <- Read10X(dataH)
#seuratTemp<-as.matrix(seuratX)[,colSums(seuratX)>500]
#write(colnames(seuratTemp),file=paste0(homeDir,"raw_files/",sample,"/CITE_fastqs/",sample,"_whitelist.txt"))

e.out<-emptyDrops(seuratX)
is.cell <- e.out$FDR <= 0.01
whitelist<-colnames(seuratX)[is.cell]
whitelist=whitelist[!is.na(whitelist)]
whitelist=gsub("-.*","",whitelist)

write(whitelist,file=paste0(homeDir,"raw_files/",sample,"/CITE_fastqs/",sample,"_whitelist.txt"))


