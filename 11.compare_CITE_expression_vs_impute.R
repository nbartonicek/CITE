library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(scImpute)
library(data.table)

homeDir="/share/ScratchGeneral/nenbar/projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4515"
sampleName="4515_CITE"
#treatment="impute"
treatment="original"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))
rObjectsDir=paste0(homeDir,"project_results/Robjects/")
system(paste0("mkdir -p ",rObjectsDir))
outDir=paste0(homeDir,"project_results/Robjects/doubletFinder/")
system(paste0("mkdir -p ",outDir))
rawDir=paste0(homeDir,"raw_files/",projectname,"/",sampleName,"/raw_data/")
imputedDir=paste0(homeDir,"raw_files/",projectname,"/",sampleName,"/imputed_data/")
system(paste0("mkdir -p ",rawDir))
system(paste0("mkdir -p ",imputedDir))

#load the annotated data
#load(paste0(homeDir,"project_results/Robjects/Human_HASHSEQ_S1_hd2_F_clean.Rdata"))
#load(paste0(rObjectsDir,"hashedCells_with_HTO.Rdata"))
#load(paste0(rObjectsDir,"imputedHashedCells_with_HTO.Rdata"))

#load("/share/ScratchGeneral/nenbar/projects/CITE/project_results/Robjects/4515_CITE_tSNE.Rdata")
load(paste0(rObjectsDir,"4515_CITE_annotated_imputed.Rdata"))


#add CITE data

cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/4515_CITE.tsv"), 
                       sep = ",", header = TRUE, row.names = 1)
rownames(cite) <- paste0("CITE_", rownames(cite))
cite <- cite[1:(dim(cite)[1]-3),]
cite <-cite[,order(colnames(cite))]
row.names(cite)[row.names(cite)=="CITE_IL-7Ralpha"]="CITE_CD127"
citeCells<-annotated_citeCells

cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)

#first make a tCell dataset
tcellL<-citeCells@ident %in% "T-cells"
imputed<-citeCells
tCellImp<-imputed@data[,tcellL]

load("/share/ScratchGeneral/nenbar/projects/CITE/project_results/Robjects/4515_CITE_tSNE.Rdata")
cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)
original<-citeCells
tCellOri<-original@data[,colnames(original@data) %in% colnames(tCellImp)]

#plot CD45 expression in ori vs imp
CD45Ori<-tCellOri["PTPRC",]
CD45Imp<-tCellImp["PTPRC",]


dfT<-data.frame(original=as.numeric(CD45Ori),imputed=as.numeric(CD45Imp))

pdf(paste0(figDir,sampleName,"_PTPRC_original_imputed_correlation.pdf"),width=6,height=6,useDingbats=FALSE)
ggplot(dfT, aes(x = original, y=imputed) ) + geom_point() +ylim(c(0,5))
dev.off()


x=as.numeric(original@assay$CITE@data[paste0("CITE_CD45"),][colnames(original@data) %in% colnames(tCellImp)])
y=as.numeric(imputed@assay$CITE@data[paste0("CITE_CD45"),][tcellL])

dfT<-data.frame(original=x,imputed=y)
     


pdf(paste0(figDir,sampleName,"_CITE_CD45_original_imputed_correlation.pdf"),width=6,height=6,useDingbats=FALSE)
ggplot(dfT, aes(x = original, y=imputed) ) + geom_point() 
dev.off()

pcaOri<-princomp(tCellOri)
pcaImp<-princomp(tCellImp)

tCellCITE=t(imputed@assay$CITE@data[,tcellL])
pcaCITE<-princomp(tCellCITE)

pdf(paste0(figDir,sampleName,"_PCA_tCells_ori.pdf"),width=5,height=5,useDingbats=FALSE)
plot(pcaOri$loadings[,1],pcaOri$loadings[,2])
dev.off()


pdf(paste0(figDir,sampleName,"_PCA_tCells_imp.pdf"),width=5,height=5,useDingbats=FALSE)
plot(pcaImp$loadings[,1],pcaImp$loadings[,2])
dev.off()

pdf(paste0(figDir,sampleName,"_PCA_tCells_CITE.pdf"),width=5,height=5,useDingbats=FALSE)
plot(pcaCITE$scores[,1],pcaCITE$scores[,2])
dev.off()














