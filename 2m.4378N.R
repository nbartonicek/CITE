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



projectname="4378N"
sampleName="4378N"
treatment="v3"

countsTotal<-list()
lossTotal<-list()
abCounts<-list()
cellNums<-list()

cat(projectname)
cat("\n")
treatment="v3"
#treatment="mt0.5_nGene7000"


matrix_dir = paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_CITE.v3/umi_count/")
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




#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))
finalObject<-paste0(RobjectsDir,"03_seurat_object_processed.RData")
#finalObject<-"/paella/TumourProgressionGroupTemp/projects/brca_mini_atlas/analysis/randomly_subsetted_20k_cells_for_testing/sDataRandomSubset.Rdata"
citeCells<-readRDS(finalObject)

cite<-cite[,colnames(cite) %in% colnames(citeCells@assays$RNA@data)]
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


citeCells[["ADT"]] <- CreateAssayObject(counts = citeAll)
citeCells <- NormalizeData(citeCells, assay = "ADT", normalization.method = "CLR")
citeCells <- ScaleData(object = citeCells, assay = "ADT")

save(citeCells,file=paste0(RobjectsDir,projectname,"_",treatment,"_CITE.Rdata"))


pdf(paste0(figDir,projectname,"_CITE_annotated_umap_garnett_seurat_cluster_call_major_PC_C_res.1.2.pdf"),width=12,height=8)
DimPlot(object = citeCells, reduction="UMAPC",label = T,group.by="garnett_seurat_cluster_call_major_PC_C_res.1.2")
dev.off()



genes<-c("nFeature_RNA","percent.mt",paste0(c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,sampleName,"_",treatment,"_feature_vst.pdf"),width=12,height=16)
FeaturePlot(object = citeCells, features = genes, ncol = 3, pt.size = 0.5,cols = c("grey", "blue"), 
            reduction = "UMAPC")

dev.off()
#}
row.names(cite)<-gsub("_","-",row.names(cite))


adt.markers <- FindAllMarkers(citeCells, assay = "ADT", only.pos = TRUE)
pdf(paste0(figDir,"CITE_garnett_seurat_cluster_call_major_PC_C_res.1.2_clusters.pdf"),width=12,height=20)
DoHeatmap(citeCells, features = unique(adt.markers$gene),raster=F, assay = "ADT", angle = 90,group.by="garnett_seurat_cluster_call_major_PC_C_res.1.2") + NoLegend()
dev.off()





pdf(paste0(figDir,"feature_plot_cite_",sampleName,"1.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[1:20], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "UMAPC", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"2.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[21:40], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "UMAPC", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"3.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[41:60], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "UMAPC", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"4.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[61:80], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "UMAPC", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"5.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[81:100], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "UMAPC", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"6.pdf"),width=12,height=5)
FeaturePlot(citeCells, features = row.names(cite)[101:111], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "UMAPC", pt.size = 0.5)
dev.off()



cite<-citeCells@assays$ADT@counts
igGs<-cite[grepl("IgG\\d",row.names(cite)),]
citeCells@meta.data$cell.size<-log1p(apply(igGs,2,sum))
citeCells@meta.data$log.median.isotype.count<-log1p(apply(igGs,2,median))

pdf(paste0(figDir,sampleName,"_",treatment,"_violin_cellSize.pdf"),width=8,height=6)
VlnPlot(object = citeCells, features="log.median.isotype.count", group.by="garnett_seurat_cluster_call_major_PC_C_res.1.2",pt.size = 0.5)
dev.off()


cite<-citeCells@assays$ADT@data
igGs<-cite[grepl("IgG\\d",row.names(cite)),]
citeCells@meta.data$median.isotype.count<-apply(igGs,2,median)

pdf(paste0(figDir,sampleName,"_",treatment,"_violin_isotype_norm.pdf"),width=8,height=6)
VlnPlot(object = citeCells, features="median.isotype.count", group.by="garnett_seurat_cluster_call_major_PC_C_res.1.2",pt.size = 0.5)
dev.off()








cite<-citeCells@assays$ADT@data
igGs<-cite[grepl("IgG\\d",row.names(cite)),]
igGsShort<-igGs[,colSums(igGs)>0]
meanIgG<-apply(igGsShort,2,mean)
cite<-citeCells@assays$ADT@counts[,colSums(igGs)>0]
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






citeCells[["ADT"]] <- CreateAssayObject(counts = citeAll)
citeCells <- NormalizeData(citeCells, assay = "ADT", normalization.method = "CLR")
citeCells <- ScaleData(object = citeCells, assay = "ADT")


adt.markers <- FindAllMarkers(citeCells, assay = "ADT", only.pos = TRUE)
pdf(paste0(figDir,"CITE_NORM_garnett_seurat_cluster_call_major_PC_C_res.1.2_clusters.pdf"),width=12,height=20)
DoHeatmap(citeCells, features = unique(adt.markers$gene),raster=F, assay = "ADT", angle = 90,group.by="garnett_seurat_cluster_call_major_PC_C_res.1.2") + NoLegend()
dev.off()





