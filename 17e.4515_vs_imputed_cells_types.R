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



projectname="4515"
sampleName="4515_CITE"
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
colnames(cite)<-paste0("CID4515_",colnames(cite))


#tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE_vst.Rdata")
#load(tSNEFile)


figDir=paste0(homeDir,"project_results/figures/impute_CITE_QC/")
system(paste0("mkdir -p ",figDir))

#annotate cells
finalObject<-paste0("/paella/TumourProgressionGroupTemp/projects/brca_mini_atlas/analysis/pathway_and_TF_analysis.DLR_working/seurat_CCA_aligned_processed.with_NewGarnett_AUCell_SCENIC.Rdata")
miniatlas<-readRDS(finalObject)

#cellannotation
annotation<-miniatlas@meta.data$int_cluster_subset_PC_C_res.1.2
annotation4515<-annotation[miniatlas@meta.data$orig.ident=="CID4515"]
barcodes<-colnames(miniatlas@assays$RNA@data)
barcodes4515<-barcodes[grepl("4515",barcodes)]
miniatlas4515<-subset(miniatlas,cells=barcodes4515)

cite <-cite[,colnames(cite) %in% colnames(miniatlas4515@assays$RNA@data) ]
miniatlas4515<-subset(miniatlas,cells=colnames(cite))

#cite <-cite[,colnames(cite) %in% colnames(citeCells@assays$RNA@data) ]
miniatlas4515[["ADT"]] <- CreateAssayObject(counts = cite)
miniatlas4515 <- NormalizeData(miniatlas4515, assay = "ADT", normalization.method = "CLR")
miniatlas4515 <- ScaleData(object = miniatlas4515, assay = "ADT")




############# add t-cell correlation ##########

adt.markers <- FindAllMarkers(miniatlas4515, assay = "ADT", only.pos = TRUE)
pdf(paste0(figDir,"4515_only_int_cluster_subset_PC_C_res.1.2_clusters.pdf"),width=8,height=6)
DoHeatmap(miniatlas4515, features = unique(adt.markers$gene),raster=F, assay = "ADT", size=3,angle = 90,group.by="int_cluster_subset_PC_C_res.1.2") + NoLegend()
dev.off()















adt.markers <- FindAllMarkers(citeCells, assay = "ADT", only.pos = TRUE)
pdf(paste0(figDir,"CITE_imputed_subsetAtlas_int_cluster_subset_PC_C_res.1.2_clusters.pdf"),width=8,height=6)
DoHeatmap(citeCells1, features = unique(adt.markers$gene),raster=F, assay = "ADT", size=4,angle = 90,group.by="int_cluster_subset_PC_C_res.1.2") + NoLegend()
dev.off()








projectname="miniatlas_all"
sampleName="miniatlas_all"
treatment="v3"

#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))
finalObject<-paste0(RobjectsDir,"03_seurat_CCA_aligned_processed.Rdata")
#finalObject<-"/paella/TumourProgressionGroupTemp/projects/brca_mini_atlas/analysis/randomly_subsetted_20k_cells_for_testing/sDataRandomSubset.Rdata"
miniatlas<-readRDS(finalObject)
citeCells1<-citeCells

cite<-cite[,colnames(cite) %in% colnames(citeCells1@assays$RNA@data)]
cite$id<-row.names(cite)
allBCs<-colnames(citeCells1@assays$RNA@data)

citeAll<-as.data.frame(matrix(0,dim(cite)[1],length(allBCs)))
row.names(citeAll)<-row.names(cite)
colnames(citeAll)<-allBCs

citeAll$id<-row.names(citeAll)
citeL<-list(cite=cite,citeAll=citeAll)

matrix.df <- ldply(citeL, melt)
sum.matrix <- acast(matrix.df, id ~ variable, sum)

citeAll<-sum.matrix
cite<-cite[,colnames(cite) %in% colnames(citeCells1@assays$RNA@data)]

citeCells1[["ADT"]] <- CreateAssayObject(counts = citeAll)
citeCells1 <- NormalizeData(citeCells1, assay = "ADT", normalization.method = "CLR")
citeCells1 <- ScaleData(object = citeCells1, assay = "ADT")


pdf(paste0(figDir,"CITE_original_miniatlas_feature_ADT.pdf"),width=12,height=16)
FeaturePlot(object = citeCells1, features = row.names(citeCells1$ADT@data), reduction="UMAPC",min.cutoff = "q05", max.cutoff = "q95", ncol = 4)
dev.off()


allCells<-colnames(citeCells1@assays$RNA@data)
nonQuery<-allCells[!(allCells %in% colnames(cite))]
reference <- subset(citeCells1,cells=colnames(cite) )
query <- subset(citeCells1,cells=nonQuery )

save(reference,file=paste0(RobjectsDir,"reference_all.Rdata"))
save(query,file=paste0(RobjectsDir,"query_all.Rdata"))


query.anchors <- FindTransferAnchors(reference = reference, query = query, project.query=T)
save(query.anchors,file=paste0(RobjectsDir,"query.anchors_all.Rdata"))


predictions1 <- TransferData(anchorset = query.anchors, refdata = reference@assays$ADT@counts)
predictionsDF<-as.data.frame(predictions1@data)
cite<-cite[order(row.names(cite)),]
predictionsDF<-predictionsDF[order(row.names(predictionsDF)),]

cite<-cbind(cite,predictionsDF)
save(cite,file=paste0(RobjectsDir,"imputed_cite_all.Rdata"))


citeCells1[["ADT"]] <- CreateAssayObject(counts = cite)
citeCells1 <- NormalizeData(citeCells1, assay = "ADT", normalization.method = "CLR")
citeCells1 <- ScaleData(object = citeCells1, assay = "ADT")

integrated.data<-citeCells1
saveRDS(integrated.data,file=paste0(RobjectsDir,"miniatlas_4515-imputed-CITE.Rdata"))

pdf(paste0(figDir,"CITE_imputed_miniatlas_feature_ADT_all.pdf"),width=12,height=16)
FeaturePlot(object = citeCells1, features = row.names(citeCells1$ADT@data), reduction="UMAPC",min.cutoff = "q05", max.cutoff = "q95", ncol = 4)
dev.off()

#citeCells1@meta.data$garnett_call_ext_major
#citeCells1 <- RenameIdents(citeCells1, citeCells1@meta.data$garnett_call_ext_major)

pdf(paste0(figDir,"CITE_annotated_umap_samples_all.pdf"),width=12,height=8)
DimPlot(object = citeCells1, reduction="UMAPC",label = T,group.by="garnett_call_ext_major")
dev.off()


pdf(paste0(figDir,"CITE_annotated_umap_samples_split_all.pdf"),width=12,height=16)
DimPlot(object = citeCells1, reduction="UMAPC",label = T,split.by="orig.ident",ncol = 4)
dev.off()



pdf(paste0(figDir,"CITE_imputed_allAtlas_garnett_call_ext_major_ridgeplot.pdf"),width=8,height=32)
RidgePlot(citeCells1, features = row.names(citeCells1$ADT@data), ncol = 2,group.by="garnett_call_ext_major")
dev.off()




pdf(paste0(figDir,"CITE_imputed_allAtlas_garnett_call_ext_major_CD3E.pdf"),width=8,height=6)
FeatureScatter(citeCells1, feature1 = "CITE-CD3", feature2 = "CD3E",group.by="garnett_call_ext_major")
dev.off()


pdf(paste0(figDir,"CITE_imputed_allAtlas_garnett_call_ext_major_CD19.pdf"),width=8,height=6)
FeatureScatter(citeCells1, feature1 = "CITE-CD19", feature2 = "CD19",group.by="garnett_call_ext_major")
dev.off()



adt.markers <- FindAllMarkers(citeCells1, assay = "ADT", only.pos = TRUE)
png(paste0(figDir,"CITE_imputed_allAtlas_garnett_call_ext_major_clusters.png"),width=8,height=6,units="in",res=75)
DoHeatmap(citeCells1, features = unique(adt.markers$gene),raster=F, assay = "ADT", angle = 90,group.by="garnett_call_ext_major") + NoLegend()
dev.off()











pdf(paste0(figDir,"CITE_annotated_umap_samples_int_cluster_subset_PC_C_res.1.2.pdf"),width=12,height=8)
DimPlot(object = citeCells1, reduction="UMAPC",label = T,group.by="int_cluster_subset_PC_C_res.1.2")
dev.off()



pdf(paste0(figDir,"CITE_imputed_allAtlas_int_cluster_subset_PC_C_res.1.2_ridgeplot.pdf"),width=9,height=40)
RidgePlot(citeCells1, features = row.names(citeCells1$ADT@data), ncol = 2,group.by="int_cluster_subset_PC_C_res.1.2")
dev.off()




pdf(paste0(figDir,"CITE_imputed_allAtlas_int_cluster_subset_PC_C_res.1.2_CD3E.pdf"),width=8,height=6)
FeatureScatter(citeCells1, feature1 = "CITE-CD3", feature2 = "CD3E",group.by="int_cluster_subset_PC_C_res.1.2")
dev.off()


pdf(paste0(figDir,"CITE_imputed_allAtlas_int_cluster_subset_PC_C_res.1.2_CD19.pdf"),width=8,height=6)
FeatureScatter(citeCells1, feature1 = "CITE-CD19", feature2 = "CD19",group.by="int_cluster_subset_PC_C_res.1.2")
dev.off()



png(paste0(figDir,"CITE_imputed_allAtlas_garnett_call_ext_major_clusters.png"),width=8,height=6,units="in",res=75)
DoHeatmap(citeCells1, features = unique(adt.markers$gene),raster=F, assay = "ADT", size=4,angle = 90,group.by="int_cluster_subset_PC_C_res.1.2") + NoLegend()
dev.off()


adt.markers <- FindAllMarkers(citeCells1, assay = "ADT", only.pos = TRUE)
png(paste0(figDir,"CITE_imputed_allAtlas_garnett_call_ext_major_clusters.png"),width=12,height=8,units="in",res=100)
DoHeatmap(citeCells1, features = unique(adt.markers$gene),raster=F, assay = "ADT", angle = 90,group.by="int_cluster_subset_PC_C_res.1.2") + NoLegend()
dev.off()


pdf(paste0(figDir,"CITE_imputed_allAtlas_int_cluster_subset_PC_C_res.1.2_ridgeplot.pdf"),width=9,height=40)
RidgePlot(tcells, features = row.names(tcells$ADT@data), ncol = 2,group.by="int_cluster_subset_PC_C_res.1.2")
dev.off()









#For 4515, TNBC and all samples calculate
#Correlation between protein and gene
#Clustering, large scale
#Overlay over clustering?





















