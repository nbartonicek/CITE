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



projectname="4497-1"
sampleName="4497-1_CITE"
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
#colnames(cite)<-paste0("CID44971_",colnames(cite))

projectname="4497-1_CHUNKS2"
sampleName="4497-1_CHUNKS2_CITE"
treatment="v2"

tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_tSNE.Rdata")
load(tSNEFile)


cite <-cite[,colnames(cite) %in% colnames(citeCells@assays$RNA@data) ]
row.names(cite)=gsub(" ","",row.names(cite))

citeCellsShort <-subset(citeCells,cells=colnames(cite))

citeCellsShort[["ADT"]] <- CreateAssayObject(counts = cite)
citeCellsShort <- NormalizeData(citeCellsShort, assay = "ADT", normalization.method = "CLR")
citeCellsShort <- ScaleData(object = citeCellsShort, assay = "ADT")


projectname="miniatlas"
sampleName="miniatlas"
treatment="v3"

#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))
#finalObject<-paste0(RobjectsDir,"03_seurat_CCA_aligned_processed.Rdata")
finalObject<-"/paella/TumourProgressionGroupTemp/projects/brca_mini_atlas/analysis/randomly_subsetted_20k_cells_for_testing/sDataRandomSubset.Rdata"
citeCells<-readRDS(finalObject)
citeCells1<-citeCells

cite<-cite[,colnames(cite) %in% colnames(citeCells1@assays$RNA@data)]
#cite$id<-row.names(cite)
allBCs<-colnames(citeCells1@assays$RNA@data)

reference <- citeCellsShort
reference <- subset(reference,features=row.names(citeCells1@assays$RNA@data))
data<-reference@assays$RNA@data
reference <- FindVariableFeatures(reference, selection.method = "vst", nfeatures = 2000)
reference <- ScaleData(reference)
reference <- RunPCA(reference, features = VariableFeatures(object = reference))
reference


query <- citeCells1
query <- subset(citeCells1,features=row.names(reference@assays$RNA@data))
#save(reference,file=paste0(RobjectsDir,"reference.Rdata"))
#save(query,file=paste0(RobjectsDir,"query.Rdata"))

query <- FindVariableFeatures(query, selection.method = "vst", nfeatures = 2000)
query <- ScaleData(query)
query <- RunPCA(query, features = VariableFeatures(object = query))
query

query.anchors <- FindTransferAnchors(reference = reference, query = query, project.query=T)
save(query.anchors,file=paste0(RobjectsDir,"4497_query.anchors.Rdata"))


predictions1 <- TransferData(anchorset = query.anchors, refdata = reference@assays$ADT@counts)
predictionsDF<-as.data.frame(predictions1@data)
cite<-cite[order(row.names(cite)),]
predictionsDF<-predictionsDF[order(row.names(predictionsDF)),]

cite<-cbind(cite,predictionsDF)
save(cite,file=paste0(RobjectsDir,"imputed_cite.Rdata"))


citeCells1[["ADT"]] <- CreateAssayObject(counts = cite)
citeCells1 <- NormalizeData(citeCells1, assay = "ADT", normalization.method = "CLR")
citeCells1 <- ScaleData(object = citeCells1, assay = "ADT")


pdf(paste0(figDir,"CITE_imputed_miniatlas_feature_ADT.pdf"),width=12,height=16)
FeaturePlot(object = citeCells1, features = row.names(citeCells1$ADT@data), reduction="UMAPC",min.cutoff = "q05", max.cutoff = "q95", ncol = 4)
dev.off()

#citeCells1@meta.data$garnett_call_ext_major
#citeCells1 <- RenameIdents(citeCells1, citeCells1@meta.data$garnett_call_ext_major)

pdf(paste0(figDir,"CITE_annotated_umap_samples.pdf"),width=12,height=8)
DimPlot(object = citeCells1, reduction="UMAPC",label = T,group.by="garnett_call_ext_major")
dev.off()


pdf(paste0(figDir,"CITE_annotated_umap_samples_split.pdf"),width=12,height=16)
DimPlot(object = citeCells1, reduction="UMAPC",label = T,split.by="orig.ident",ncol = 4)
dev.off()



pdf(paste0(figDir,"CITE_annotated_umap_samples_4515_CITE.pdf"),width=12,height=16)
FeaturePlot(object = citeCells1, reduction="UMAPC",features = row.names(citeCells1$ADT@data),label = T,split.by="orig.ident",ncol = 4)
dev.off()


#
pdf(paste0(figDir,"CITE_imputed_subsetAtlas_garnett_call_ext_major_ridgeplot.pdf"),width=8,height=32)

RidgePlot(citeCells1, features = row.names(citeCells1$ADT@data), ncol = 2,group.by="garnett_call_ext_major")

dev.off()


pdf(paste0(figDir,"CITE_imputed_subsetAtlas_garnett_call_ext_major_CD3E.pdf"),width=8,height=6)
FeatureScatter(citeCells1, feature1 = "CITE-CD3", feature2 = "CD3E",group.by="garnett_call_ext_major")
dev.off()


pdf(paste0(figDir,"CITE_imputed_subsetAtlas_garnett_call_ext_major_CD19.pdf"),width=8,height=6)
FeatureScatter(citeCells1, feature1 = "CITE-CD19", feature2 = "CD19",group.by="garnett_call_ext_major")
dev.off()

adt.markers <- FindAllMarkers(citeCells1, assay = "ADT", only.pos = TRUE)
pdf(paste0(figDir,"CITE_imputed_subsetAtlas_garnett_call_ext_major_clusters.pdf"),width=8,height=6)
DoHeatmap(citeCells1, features = unique(adt.markers$gene),raster=F, assay = "ADT", angle = 90,group.by="garnett_call_ext_major") + NoLegend()
dev.off()





pdf(paste0(figDir,"CITE_annotated_umap_samples_int_cluster_subset_PC_C_res.1.2.pdf"),width=12,height=8)
DimPlot(object = citeCells1, reduction="UMAPC",label = T,group.by="int_cluster_subset_PC_C_res.1.2")
dev.off()



pdf(paste0(figDir,"CITE_imputed_subsetAtlas_int_cluster_subset_PC_C_res.1.2_ridgeplot.pdf"),width=8,height=40)
RidgePlot(citeCells1, features = row.names(citeCells1$ADT@data), ncol = 2,group.by="int_cluster_subset_PC_C_res.1.2")
dev.off()




pdf(paste0(figDir,"CITE_imputed_subsetAtlas_int_cluster_subset_PC_C_res.1.2_CD3E.pdf"),width=8,height=6)
FeatureScatter(citeCells1, feature1 = "CITE-CD3", feature2 = "CD3E",group.by="int_cluster_subset_PC_C_res.1.2")
dev.off()


pdf(paste0(figDir,"CITE_imputed_subsetAtlas_int_cluster_subset_PC_C_res.1.2_CD19.pdf"),width=8,height=6)
FeatureScatter(citeCells1, feature1 = "CITE-CD19", feature2 = "CD19",group.by="int_cluster_subset_PC_C_res.1.2")
dev.off()



adt.markers <- FindAllMarkers(citeCells1, assay = "ADT", only.pos = TRUE)
pdf(paste0(figDir,"CITE_imputed_subsetAtlas_int_cluster_subset_PC_C_res.1.2_clusters.pdf"),width=8,height=6)
DoHeatmap(citeCells1, features = unique(adt.markers$gene),raster=F, assay = "ADT", size=4,angle = 90,group.by="int_cluster_subset_PC_C_res.1.2") + NoLegend()
dev.off()














