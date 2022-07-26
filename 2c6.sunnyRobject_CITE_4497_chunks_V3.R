library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(plyr)
library(reshape2)

#library(DropletUtils)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4497-1_CHUNKS2"
sampleName="4497-1_CHUNKS2_CITE"
treatment="v3"
#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))
rawDir=paste0(homeDir,"raw_files/",projectname,"/",sampleName,"/raw_data/")
system(paste0("mkdir -p ",rawDir))

cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,".tsv"), 
                       sep = ",", header = TRUE, row.names = 1)
rownames(cite) <- paste0("CITE_", rownames(cite))
cite <- cite[1:(dim(cite)[1]-3),]
cite <-cite[,order(colnames(cite))]
colnames(cite)<-paste0("CID44971CHUNKS2_",colnames(cite))

seuratFile="/share/ScratchGeneral/sunwu/projects/cryopreservation_paper/seurat_v3/07_CITESeq/01_annotation_CID44971CHUNKS2/02_Rdata/04_seurat_object_processed_filtered_annotated.RData"
citeCells<-readRDS(seuratFile)

#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells@assays$RNA@data) ]

#adjust for the missing cells - make them 0
cite$id<-row.names(cite)
allBCs<-colnames(citeCells@assays$RNA@data)

citeAll<-as.data.frame(matrix(0,dim(cite)[1],length(allBCs)))
row.names(citeAll)<-row.names(cite)
colnames(citeAll)<-allBCs

citeAll$id<-row.names(citeAll)
citeL<-list(cite=cite,citeAll=citeAll)

matrix.df <- ldply(citeL, melt)
sum.matrix <- acast(matrix.df, id ~ variable, sum)

citeAll<-sum.matrix
cite<-citeAll

#add data to the main object
row.names(cite)<-gsub("_","-",row.names(cite))
citeCells[["ADT"]] <- CreateAssayObject(counts = cite)
citeCells <- NormalizeData(citeCells, assay = "ADT", normalization.method = "CLR")
citeCells <- ScaleData(object = citeCells, assay = "ADT")

#find markers that work using clustering
adt.markers <- FindAllMarkers(citeCells, assay = "ADT", only.pos = TRUE)

#on them do the ridge, heatmap and feature
short.markers.abs<-c("CITE-CD2","CITE-CTLA-4","CITE-MCAM","CITE-CD34","CITE-CD31","CITE-CD4","CITE-CD8a","CITE-CD11C","CITE-CD15")

pdf(paste0(figDir,"ridgeplot_forSunny_markers_",sampleName,".pdf"),width=12,height=30)

RidgePlot(citeCells, features = adt.markers$gene, ncol = 3)

dev.off()



pdf(paste0(figDir,"ridgeplot_forSunny_",sampleName,".pdf"),width=12,height=12)

RidgePlot(citeCells, features = short.markers.abs, ncol = 3)

dev.off()



pdf(paste0(figDir,"CITE_clusters_",sampleName,".pdf"),width=12,height=20)
p<-DoHeatmap(citeCells, features = unique(adt.markers$gene),raster=F, assay = "ADT", angle = 90) + NoLegend()
print(p)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"1.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[1:20], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"ridgeplot_",sampleName,".pdf"),width=12,height=4)

RidgePlot(citeCells, features = row.names(cite), cols = cols,ncol = 3)

dev.off()





