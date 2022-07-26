library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(sctransform)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(plyr)

homeDir="../../../projects/CITE/"
RobjectsDir=paste0(homeDir,"project_results/Robjects/")



projectname="4378N_citenorm"
sampleName="4378N_citenorm"
treatment="v3"

countsTotal<-list()
lossTotal<-list()
abCounts<-list()
cellNums<-list()

cat(projectname)
cat("\n")
treatment="v3"
#treatment="mt0.5_nGene7000"


matrix_dir = paste0(homeDir,"raw_files/4378N/CITE_fastqs/4378N_CITE.v3/umi_count/")
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
cite<-as.data.frame(mat)


rownames(cite) <- paste0("CITE_", rownames(cite))
rownames(cite) <- gsub("-\\w{15}$","",rownames(cite))
cite <-cite[,order(colnames(cite))]
cite <-cite[!grepl("unmapped",row.names(cite)), ]
cite <-cite[rowSums(cite)>1000,]
colnames(cite)<-paste0("CID4378N_",colnames(cite))




projectname="miniatlas_normal_citenorm"
sampleName="miniatlas_normal_citenorm"
treatment="v3"

#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))
finalObject<-paste0(RobjectsDir,"03_seurat_CCA_aligned_processed_normal.Rdata")
#finalObject<-"/paella/TumourProgressionGroupTemp/projects/brca_mini_atlas/analysis/randomly_subsetted_20k_cells_for_testing/sDataRandomSubset.Rdata"
citeCells<-readRDS(finalObject)
citeCells1<-citeCells

cite<-cite[,colnames(cite) %in% colnames(citeCells1@assays$RNA@data)]

igGs<-cite[grepl("IgG\\d",row.names(cite)),]
igGsShort<-igGs[,colSums(igGs)>0]
meanIgG<-apply(igGsShort,2,mean)
cite<-cite[,colSums(igGs)>0]
cite<-as.data.frame(cite)

citeAllNorm<-t(apply(cite,1,function(x){as.integer(x/meanIgG)}))
colnames(citeAllNorm)<-colnames(cite)
citeAllNorm<-as.data.frame(citeAllNorm)


citeAllNorm$id<-row.names(citeAllNorm)
allBCs<-colnames(citeCells@assays$RNA@data)

citeAll<-as.data.frame(matrix(0,dim(citeAllNorm)[1],length(allBCs)))
row.names(citeAll)<-row.names(citeAllNorm)
colnames(citeAll)<-allBCs

citeAll$id<-row.names(citeAll)
citeL<-list(cite=citeAllNorm,citeAll=citeAll)

matrix.df <- ldply(citeL, melt)
sum.matrix <- acast(matrix.df, id ~ variable, sum)

citeAll<-sum.matrix



citeCells1[["ADT"]] <- CreateAssayObject(counts = citeAll)
citeCells1 <- NormalizeData(citeCells1, assay = "ADT", normalization.method = "CLR")
citeCells1 <- ScaleData(object = citeCells1, assay = "ADT")


allCells<-colnames(citeCells1@assays$RNA@data)
nonQuery<-allCells[!(allCells %in% colnames(cite))]
reference <- subset(citeCells1,cells=colnames(cite) )
query <- subset(citeCells1,cells=nonQuery )

#save(reference,file=paste0(RobjectsDir,"reference_normal.Rdata"))
#save(query,file=paste0(RobjectsDir,"query_normal.Rdata"))


query.anchors <- FindTransferAnchors(reference = reference, query = query, project.query=T)
save(query.anchors,file=paste0(RobjectsDir,"query.anchors_normal.Rdata"))


predictions1 <- TransferData(anchorset = query.anchors, refdata = reference@assays$ADT@counts)
predictionsDF<-as.data.frame(predictions1@data)
cite<-cite[order(row.names(cite)),]
predictionsDF<-predictionsDF[order(row.names(predictionsDF)),]

cite<-cbind(cite,predictionsDF)
save(cite,file=paste0(RobjectsDir,"imputed_cite_normal_citenormal.Rdata"))


citeCells1[["ADT"]] <- CreateAssayObject(counts = cite)
citeCells1 <- NormalizeData(citeCells1, assay = "ADT", normalization.method = "CLR")
citeCells1 <- ScaleData(object = citeCells1, assay = "ADT")

integrated.data<-citeCells1
saveRDS(integrated.data,file=paste0(RobjectsDir,"miniatlas_normal_4378N-imputed-CITE_citenormal.Rdata"))

row.names(cite)<-gsub("_","-",row.names(cite))

pdf(paste0(figDir,"feature_plot_cite_",sampleName,"1.pdf"),width=12,height=10)
FeaturePlot(citeCells1, features = row.names(cite)[1:20], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "UMAPC", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"2.pdf"),width=12,height=10)
FeaturePlot(citeCells1, features = row.names(cite)[21:40], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "UMAPC", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"3.pdf"),width=12,height=10)
FeaturePlot(citeCells1, features = row.names(cite)[41:60], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "UMAPC", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"4.pdf"),width=12,height=10)
FeaturePlot(citeCells1, features = row.names(cite)[61:80], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "UMAPC", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"5.pdf"),width=12,height=10)
FeaturePlot(citeCells1, features = row.names(cite)[81:100], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "UMAPC", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"6.pdf"),width=12,height=5)
FeaturePlot(citeCells1, features = row.names(cite)[101:111], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "UMAPC", pt.size = 0.5)
dev.off()


adt.markers <- FindAllMarkers(citeCells1, assay = "ADT", only.pos = TRUE)
pdf(paste0(figDir,"CITE_garnett_call_ext_major_clusters.pdf"),width=12,height=20)
DoHeatmap(citeCells1, features = unique(adt.markers$gene),raster=F, assay = "ADT", angle = 90,group.by="garnett_call_ext_major") + NoLegend()
dev.off()


pdf(paste0(figDir,"CITE_imputed_miniatlas_normal_garnett_call_ext_major_ridgeplot_citenormal.pdf"),width=8,height=32)
RidgePlot(citeCells1, features = row.names(citeCells1$ADT@data), ncol = 2,group.by="garnett_call_ext_major")
dev.off()















