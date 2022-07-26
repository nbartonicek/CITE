library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(R.utils)
library(scImpute)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(reshape2)

dims=40
resolution=0.8
perplexities=c(5,10,15,20,50)
perplexity=20


conditions<-paste(dims,resolution,perplexity,sep="_")
cat(conditions)

homeDir="/share/ScratchGeneral/nenbar/projects/CITE/"
homeDir="../../../projects/CITE/"

projectname="4515"
sampleName="4515_CITE"
treatment="impute"

rawDir=paste0(homeDir,"raw_files/",projectname,"/",sampleName,"/raw_data/")
imputedDir=paste0(homeDir,"raw_files/",projectname,"/",sampleName,"/imputed_data/")
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_gating/")
system(paste0("mkdir -p ",rawDir))
system(paste0("mkdir -p ",imputedDir))
system(paste0("mkdir -p ",figDir))

cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/4515_CITE.tsv"), 
                       sep = ",", header = TRUE, row.names = 1)
rownames(cite) <- paste0("CITE_", rownames(cite))
cite <- cite[1:(dim(cite)[1]-3),]
cite <-cite[,order(colnames(cite))]
row.names(cite)[row.names(cite)=="CITE_IL-7Ralpha"]="CITE_CD127"

#normalize by IgG
igGs<-read.table(paste0(homeDir,"project_results/Robjects/igG_content.csv"),sep=",",header=F)
temp=cite
#normalize each CITE by corresponding IgG
for(i in 1:dim(igGs)[1]){

  ab=as.character(igGs[i,1])
  cat(ab)
  cat("\n")
  iGg<-as.character(igGs[i,2])
  cite[ab,]=pmax(0,c(as.numeric(cite[ab,])-as.numeric(cite[iGg,])))

}



tSNEFile=paste0(imputedDir,sampleName,"_tSNE.Rdata")

load(tSNEFile)



n <- length(unique(citeCells@ident))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F",cols)

pdf(paste0(figDir,sampleName,"_tsne.pdf"),width=12,height=8)
TSNEPlot(object = citeCells, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         colors.use = cols,
         no.axes=T
)

dev.off()



current.cluster.ids <- c(sort(unique(citeCells@ident)))-1
new.cluster.ids <- c(
  "epithelial_1", 
  "unknown epithelial\nCD15+_1", 
  "B-cells\nmemory", 
  "epithelial_2",
  "epithelial_3", 
  "dendritic-macrophage\nmacrophage-monocyte", 
  "T-cells", 
  "B-cells\nnaive",
  "APC",
  "NK",
  "CAF",
  "endothelial",
  "plasmablast",
  "pericytes",
  "unknown epithelial\nCD15+_2",
  "unknown epithelial",
  "epithelial_4",
  "unknown\nlow gene",
  "epithelial_5",
  "tumour_8",
  "unknown"
  )
annotated_citeCells=citeCells
annotated_citeCells@ident <- plyr::mapvalues(x = annotated_citeCells@ident, from = current.cluster.ids, to = new.cluster.ids)
save(annotated_citeCells,file=paste0(rObjectsDir,sampleName,"_annotated_imputed.Rdata"))

pdf(paste0(figDir,sampleName,"_annotated_tsne.pdf"),width=12,height=8,useDingbats=FALSE)
TSNEPlot(object = annotated_citeCells, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         colors.use = cols,
         no.axes=T
)

dev.off()



#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)


pdf(paste0(figDir,"feature_plot_Norm_cite_",sampleName,".pdf"),width=16,height=10,useDingbats=FALSE)
FeaturePlot(citeCells, features.plot = row.names(cite), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 6, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()

#Featureplot expression + CITE

selectedCITEs<-c("CD3E","CD19","CD14","FCGR3A","CD86","PTPRC",
  paste0("CITE_",c("CD3","CD19","CD14","CD16","CD86","CD45")))
pdf(paste0(figDir,"feature_plot_tx_vs_CITE_",sampleName,".pdf"),width=14,height=4,useDingbats=FALSE)
FeaturePlot(citeCells, features.plot = selectedCITEs, min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 6, cols.use = c("lightgrey", "blue"), pt.size = 0.5,no.axes=T)
dev.off()


selectedCITEs<-c("CD274","FUT4",
  paste0("CITE_",c("PD-L1","CD15")))
pdf(paste0(figDir,"feature_plot_tx_vs_CITE_addon_",sampleName,".pdf"),width=4,height=4,useDingbats=FALSE)
FeaturePlot(citeCells, features.plot = selectedCITEs, min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 2, cols.use = c("lightgrey", "blue"), pt.size = 0.5,no.axes=T)
dev.off()

pdf(paste0(testFigDir,"ridgeplot.pdf"),width=12,height=30)
RidgePlot(citeCells, features.plot = selectedCITEs, cols.use = cols,nCol = 4)
dev.off()


genes<-c("nGene","percent.mito","ACTB","EPCAM","KRT8","KRT5","COL1A1","PDGFRA","PDGFRA","ACTA2","CD34","PECAM1","PTPRC","CD68","FCGR3A","CD19","CD3D","FOXP3")
pdf(paste0(figDir,sampleName,"_feature.pdf"),width=16,height=8,useDingbats=FALSE)
FeaturePlot(object = citeCells, features.plot = genes, nCol = 6, pt.size = 0.5,cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

dev.off()


genes<-c("nGene","percent.mito","HLA-DQA1","IL2RA","IGHD","IGHG1","TGFB1","THBD","ERBB2","ESR1","AR","ITGAE","ITGAX","PRF1","GZMB","GZMA","NKG7","GNLY","TRGC1","IGHM","IGHA1","TRAC","TRBC1","CD1C","IFNG","IL3RA","CLEC4C","IL10","CD163","CXCR3","TRDC","CSF1R","STC1","CD38","TRGC2","HLA-A","CD4","CD8A","CD274")
pdf(paste0(figDir,sampleName,"_feature2.pdf"),width=16,height=14,useDingbats=FALSE)
FeaturePlot(object = citeCells, features.plot = genes, nCol = 6, pt.size = 0.5,cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

dev.off()


genes<-c("nGene","percent.mito","ITGAM","ITGAE","IFNG","CXCL10","BATF3","CXCL9","ADGRE1")
pdf(paste0(figDir,sampleName,"_feature3.pdf"),width=8,height=8,useDingbats=FALSE)
FeaturePlot(object = citeCells, features.plot = genes, nCol = 3, pt.size = 0.5,cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

dev.off()



immune=c(2,5,6,7,8,9,12)
tcell=c(6)
myeloid=c(5,8)

immuneL<-citeCells@ident %in% immune
tcellL<-citeCells@ident %in% tcell
myeloidL<-citeCells@ident %in% myeloid


types<-list()
types[["immune"]]<-immuneL
types[["T-cell"]]<-tcellL
types[["myeloid"]]<-myeloidL


combination=list()
combination[["immune"]]=list()
combination[["immune"]][[1]]=c("CD3","CD19")
combination[["T-cell"]]=list()
combination[["T-cell"]][[1]]=c("CD4","CD8a")
combination[["T-cell"]][[2]]=c("CD45RA","CD45RO")
combination[["T-cell"]][[3]]=c("CD4","CD25")
combination[["T-cell"]][[4]]=c("CD127","CTLA-4")
combination[["T-cell"]][[5]]=c("PD-1","CD28")
combination[["myeloid"]]=list()
combination[["myeloid"]][[1]]=c("CD86","CD40")

for(type in names(types)){

  for(i in 1:length(combination[[type]])){
      comb=combination[[type]][[i]]
      x=as.numeric(citeCells@assay$CITE@data[paste0("CITE_",comb[1]),][types[[type]]])
      y=as.numeric(citeCells@assay$CITE@data[paste0("CITE_",comb[2]),][types[[type]]])
      cat(length(x))
      cat("\n")
      x=jitter(x,500)
      y=jitter(y,500)
      pdf(paste0(figDir,"GenePlot_",comb[1],"_",comb[2],".pdf"),width=8,height=8)
      #GenePlot(citeCells,paste0("CITE_",ab1),paste0("CITE_",ab2),use.raw=F,use.scaled=T,cex.use = 0.6,col.use = "black")
      plot(x,y,pch=16,xlab=comb[1],ylab=comb[2],cex.lab=1.4,main=type,col="red",cex=1.2)
      dev.off()

  }

}
resolution=1.4

pcsData<-as.data.frame(citeCells@dr$pca@jackstraw@overall.p.values)
pcs<-pcsData[pcsData$Score<0.01,"PC"]

citeCellsRes <- FindClusters(citeCells, dims.use = pcs, resolution = resolution,print.output = FALSE,force.recalc=T)
citeCellsRes <- RunTSNE(citeCellsRes, dims.use = pcs, perplexity=perplexity)


pdf(paste0(figDir,sampleName,"_Res1.4_tsne.pdf"),width=12,height=8)
TSNEPlot(object = citeCellsRes, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         colors.use = cols
)
dev.off()

cluster10.markers <- FindMarkers(object = citeCellsRes, ident.1 = 10, ident.2 = 11, min.pct = 0.25)
print(x = head(x = cluster10.markers, n = 30))
cluster11.markers <- FindMarkers(object = citeCellsRes, ident.1 = 11, ident.2 = 10, min.pct = 0.25)
print(x = head(x = cluster11.markers, n = 30))
cluster11.markers <- FindMarkers(object = citeCellsRes, ident.1 = 8, ident.2 = c(10,11), min.pct = 0.25)
print(x = head(x = cluster11.markers, n = 30))

pdf(paste0(figDir,sampleName,"_violin.pdf"),width=16,height=6,useDingbats=FALSE)

VlnPlot(object = citeCells, ident.include=6,features.plot = c("CITE_CD28","CITE_B7-H4","CITE_CD40","CITE_CD86","CITE_PD-1","CITE_PD-L1","CITE_CTLA-4","CITE_TIGIT"), use.raw = F, nCol=8)

dev.off()

pdf(paste0(figDir,sampleName,"_ridgeCD3.pdf"),width=16,height=6,useDingbats=FALSE)
RidgePlot(citeCells, features.plot = c("CD3G","CITE_CD3"), ident.include=6,cols.use = cols,nCol = 2,y.max=2.5)
dev.off()

cellTypes<-as.character(annotated_citeCells@ident)

cellTypes[grepl("epithelial",cellTypes)]="epithelial"
cellTypes[grepl("tumour",cellTypes)]="epithelial"

cellTypes[grepl("B-cells",cellTypes)]="B-cells"
cellTypes[grepl("epithelial",cellTypes)]="epithelial"
cellTypes[grepl("APC",cellTypes)]="myeloid"
cellTypes[grepl("macrophage",cellTypes)]="myeloid"
cellTypes[grepl("unknown",cellTypes)]="unknown"


df=data.frame(cellType=cellTypes,id=1:length(cellTypes))
colsR=rev(cols)
sortedT<-sort(table(cellTypes),decreasing=F)

dfT<-data.frame(counts=as.numeric(sortedT),cellType=names(sortedT))
dfT$cellType=factor(dfT$cellType,levels=names(sortedT))

pdf(paste0(figDir,sampleName,"_pie.pdf"),width=16,height=6,useDingbats=FALSE)
ggplot(dfT, aes(x = factor(1), y=counts,fill=factor(cellType)) ) + geom_bar(width = 1,stat="identity")+coord_polar(theta = "y") + scale_fill_manual(values = colsR)+ guides(col = guide_legend(reverse = TRUE))
#ggplot(df,aes(cellType))+geom_bar()+ coord_polar(theta="y")
dev.off()


pdf(paste0(figDir,sampleName,"_GenePlot_CD45_nonimputed.pdf"),width=16,height=6,useDingbats=FALSE)
GenePlot(citeCells,"PTPRC","CITE_CD45",cell.ids=colnames(citeCells@data)[tcellL],col="black")
dev.off()

#load("/share/ScratchGeneral/nenbar/projects/CITE/project_results/Robjects/4515_CITE_tSNE.Rdata")
