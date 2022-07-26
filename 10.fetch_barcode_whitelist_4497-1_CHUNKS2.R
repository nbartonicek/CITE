library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(DropletUtils)


homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4497-1_CHUNKS2"
sample="4497-1_CHUNKS2"
treatment="tumor"
figDir=paste0(homeDir,"project_results/figures/",sample,"/")
system(paste0("mkdir -p ",figDir))
rObjectsDir=paste0(homeDir,"project_results/Robjects/")
system(paste0("mkdir -p ",figDir))

#emptyDrops
sample="4497-2_CELLSUS"

dataH <- paste0(homeDir,"raw_files/",sample,"/outs/raw_gene_bc_matrices/GRCh38")
dataM <- paste0(homeDir,"raw_files/",sample,"/outs/raw_gene_bc_matrices/mm10")

seuratH <- Read10X(dataH)
seuratM <- Read10X(dataM)
seuratX<-rbind(seuratH,seuratM)
seuratTemp<-as.matrix(seuratX[,colSums(seuratX)>500])
write(colnames(seuratTemp),file=paste0(homeDir,"raw_files/",sample,"/CITE_fastqs/",sample,"_whitelist.txt"))

#this produces 3262 cells (cutoff 500)
e.out<-emptyDrops(seuratX)
is.cell <- e.out$FDR <= 0.01
whitelist<-colnames(seuratX)[is.cell]
whitelist=whitelist[!is.na(whitelist)]
whitelist=gsub("-.*","",whitelist)

write(whitelist,file=paste0(homeDir,"raw_files/",sample,"/CITE_fastqs/",sample,"_whitelist.txt"))


#emptyDrops
sample="4497-1_CHUNKS2"

dataH <- paste0(homeDir,"raw_files/",sample,"/outs/raw_gene_bc_matrices/GRCh38")
dataM <- paste0(homeDir,"raw_files/",sample,"/outs/raw_gene_bc_matrices/mm10")

seuratH <- Read10X(dataH)
seuratM <- Read10X(dataM)
seuratX<-rbind(seuratH,seuratM)
seuratX<-rbind(seuratH,seuratM)
seuratTemp<-as.matrix(seuratX[,colSums(seuratX)>500])
write(colnames(seuratTemp),file=paste0(homeDir,"raw_files/",sample,"/CITE_fastqs/",sample,"_whitelist.txt"))

e.out<-emptyDrops(seuratX)
is.cell <- e.out$FDR <= 0.01
whitelist<-colnames(seuratX)[is.cell]
whitelist=whitelist[!is.na(whitelist)]
whitelist=gsub("-.*","",whitelist)

write(whitelist,file=paste0(homeDir,"raw_files/",sample,"/CITE_fastqs/",sample,"_whitelist.txt"))


