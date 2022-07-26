library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(sctransform)
library(ggplot2)
library(descend)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4515"
sampleName="4515_CITE"
treatment="v3"
#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))


matrix_dir = paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/4515_CITE.v3.tsv/umi_count/")
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "features.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
cite<-as.data.frame(mat)

rownames(cite) <- paste0("CITE_", rownames(cite))
rownames(cite) <- gsub("-\\w{15}$","",rownames(cite))
cite <-cite[,order(colnames(cite))]

tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE_vst.Rdata")

#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)

