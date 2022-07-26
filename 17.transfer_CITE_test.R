library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(sctransform)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(cowplot)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/QC/")




projectname="4497-2"
sampleName="4497-2_CELLSUS"
treatment="v3"

countsTotal<-list()
lossTotal<-list()
abCounts<-list()
cellNums<-list()

cat(projectname)
cat("\n")
treatment="v3"
#treatment="mt0.5_nGene7000"
system(paste0("mkdir -p ",figDir))
citeDir=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,".v3")

finalObject<-paste0("../project_results/Robjects/",sampleName,"_",treatment,"_CITE.Rdata")
load(finalObject)
citeCells1<-citeCells

#eliminate mouse and low quality
citeCells1Counts<-citeCells1@assays$RNA@counts
humanIdx<-grepl("GRCh38-",row.names(citeCells1Counts))
speciesRatio<-apply(citeCells1Counts,2,function(x){sum(x[humanIdx])/sum(x)})
citeCells1Counts<-citeCells1Counts[,speciesRatio>=0.9]

colnames(citeCells1Counts)=paste0(colnames(citeCells1Counts),"_1")


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

colnames(cite)<-paste0(colnames(cite),"_1")
cite1<-cite[colnames(cite) %in% colnames(citeCells1Counts)]

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
system(paste0("mkdir -p ",figDir))
citeDir=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,".v3")

finalObject<-paste0("../project_results/Robjects/",sampleName,"_",treatment,"_CITE.Rdata")
load(finalObject)

citeCells2<-citeCells
citeCells2Counts<-citeCells2@assays$RNA@counts
colnames(citeCells2Counts)=paste0(colnames(citeCells2Counts),"_2")


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
colnames(cite)<-paste0(colnames(cite),"_2")
cite2<-cite[colnames(cite) %in% colnames(citeCells2Counts)]

######## integrate ########
citeCells1Counts<-citeCells1Counts[grepl("GRCh38",row.names(citeCells1Counts)),]
row.names(citeCells1Counts)<-gsub("GRCh38-","",row.names(citeCells1Counts))
citeCells1Counts<-citeCells1Counts[row.names(citeCells1Counts) %in% row.names(citeCells2Counts),]
citeCells2Counts<-citeCells2Counts[row.names(citeCells2Counts) %in% row.names(citeCells1Counts),]

c1<-CreateSeuratObject(citeCells1Counts)
c2<-CreateSeuratObject(citeCells2Counts)
integrated<-list()

integrated[["4497"]]<-NormalizeData(c1,verbose=FALSE)
integrated[["4515"]]<-NormalizeData(c2,verbose=FALSE)
integrated[["4497"]] <- FindVariableFeatures(integrated[["4497"]], selection.method = "vst", nfeatures = 2000, 
        verbose = FALSE)
integrated[["4515"]] <- FindVariableFeatures(integrated[["4515"]], selection.method = "vst", nfeatures = 2000, 
        verbose = FALSE)
#integrated[["4497"]][["ADT"]] <- CreateAssayObject(counts = cite1)
#integrated[["4497"]] <- NormalizeData(integrated[["4497"]], assay = "ADT", normalization.method = "CLR")
#integrated[["4497"]] <- ScaleData(object = integrated[["4497"]], assay = "ADT")

integrated[["4515"]][["ADT"]] <- CreateAssayObject(counts = cite2)
integrated[["4515"]] <- NormalizeData(integrated[["4515"]], assay = "ADT", normalization.method = "CLR")
integrated[["4515"]] <- ScaleData(object = integrated[["4515"]], assay = "ADT")

save(integrated,file=paste0("../project_results/Robjects/CITE_integrated_transfer_test.Rdata"))

anchors <- FindIntegrationAnchors(object.list = integrated, dims = 1:30)
cite.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)


DefaultAssay(cite.integrated) <- "integrated"


cite.integrated <- ScaleData(cite.integrated, verbose = FALSE)
cite.integrated <- RunPCA(cite.integrated, npcs = 30, verbose = FALSE)
cite.integrated <- RunUMAP(cite.integrated, reduction = "pca", dims = 1:30)
cite.integrated <- FindNeighbors(cite.integrated, dims.use = 1:30, resolution = 0.8)
cite.integrated <- FindClusters(cite.integrated, dims.use = 1:30, resolution = 0.8)
cite.integrated <- RunTSNE(cite.integrated, dims.use = 1:30,check_duplicates=FALSE)

cite.integrated@meta.data$orig.ident=rep(c("4497","4515"),times=c(dim(citeCells1Counts)[2],dim(citeCells2Counts)[2]))

save(cite.integrated,file=paste0("../project_results/Robjects/CITE_transfer_test.Rdata"))

pdf(paste0(figDir,"CITE_transfer_test_bothsamples.pdf"),width=12,height=8)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F","darkblue","orange")
p1 <- DimPlot(cite.integrated, reduction = "tsne", group.by = "orig.ident",pt.size = 1.5,cols=cols)
p1
dev.off()


query <- integrated[["4497"]]
query.anchors <- FindTransferAnchors(reference = cite.integrated, query = query, 
    dims = 1:30)
predictions1 <- TransferData(anchorset = query.anchors, refdata = cite.integrated@assays$ADT@counts, 
    dims = 1:30)

cites<-row.names(predictions1@data)

pdf(paste0(figDir,"imputed_4497.pdf"),width=12,height=12)
  par(mfrow=c(5,5))
  for(ab in cites){
    cat(ab)
    cat("\n")
    df4515=log1p(as.numeric(predictions1@data[ab,]))
    abUnd<-gsub("CITE-","CITE_",ab)
    df4497=log1p(as.numeric(cite1[abUnd,]))
    corEstimate<-as.numeric(cor.test(df4515,df4497)$estimate)
    corEstimate<-sprintf("%1.5f",corEstimate)
    plot(df4515,df4497,main=paste0(ab," cor: ",corEstimate))

  }

dev.off()


imputed4497<-integrated[["4497"]]
imputed4497[["ADT"]] <- CreateAssayObject(counts = cite1)
imputed4497 <- NormalizeData(imputed4497, assay = "ADT", normalization.method = "CLR")
imputed4497 <- ScaleData(object = imputed4497, assay = "ADT")
clr4497<-as.data.frame(imputed4497@assays$ADT@data)


imputed4497<-integrated[["4497"]]
imputed4497[["ADT"]] <- CreateAssayObject(counts = predictions1@data)
imputed4497 <- NormalizeData(imputed4497, assay = "ADT", normalization.method = "CLR")
imputed4497 <- ScaleData(object = imputed4497, assay = "ADT")
imputedCITE<-as.data.frame(imputed4497@assays$ADT@data)

cites<-row.names(predictions1@data)


pdf(paste0(figDir,"imputed_CLR_4497.pdf"),width=12,height=12)
  par(mfrow=c(5,5))
  for(ab in cites){
    cat(ab)
    cat("\n")
    df4515=as.numeric(imputedCITE[ab,])
    df4497=as.numeric(clr4497[ab,])
    corEstimate<-as.numeric(cor.test(df4515,df4497)$estimate)
    corEstimate<-sprintf("%1.5f",corEstimate)
    plot(df4515,df4497,main=paste0(ab," cor: ",corEstimate))

  }

dev.off()



pdf(paste0(figDir,"CITE_integrated_tsne_clusters.pdf"),width=12,height=8)

p1 <- DimPlot(cite.integrated, reduction = "tsne")
p1
dev.off()




pdf(paste0(figDir,"CITE_integrated_tsne.pdf"),width=12,height=8)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F","darkblue","orange")
DimPlot(object = cite.integrated, 
         label = T, 
         pt.size = 2, 
         label.size = 6,
         cols = cols
)
dev.off()


genes<-c("nFeature_RNA","percent.mt",paste0(c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,"CITE_integrated_feature_genes.pdf"),width=12,height=16)
FeaturePlot(object = cite.integrated, features = genes, ncol = 3, pt.size = 0.5,cols = c("grey", "blue"), 
            reduction = "tsne")

dev.off()
#}


pdf(paste0(figDir,"CITE_integrated_feature_ADT.pdf"),width=12,height=16)
FeaturePlot(object = cite.integrated, features = row.names(cite.integrated$ADT@data), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)
dev.off()

############ plot the correlations only for T cells #################


bcs<-c(1:dim(predictions@data)[2])[cite.integrated@meta.data$seurat_clusters %in% c(0,1,3,6)]
bcs<-bcs[!is.na(bcs)]

cites<-row.names(predictions@data)
#cites<-cites[!grepl("IL-7",cites)]
#cites<-cites[!grepl("CTLA-4",cites)]
#cites<-cites[!grepl("PD-L1",cites)]
#cites<-cites[!grepl("PD-1",cites)]
#cites<-cites[!grepl("B7-H4",cites)]

pdf(paste0(figDir,"imputed_Tcells_4497.pdf"),width=12,height=12)
  par(mfrow=c(5,5))
  for(ab in cites){
    cat(ab)
    cat("\n")
    df4515=as.numeric(imputedCITE[ab,])[bcs]
    df4497=as.numeric(clr4497[ab,])[bcs]
    corEstimate<-as.numeric(cor.test(df4515,df4497)$estimate)
    corEstimate<-sprintf("%1.5f",corEstimate)
    plot(df4515,df4497,main=paste0(ab," cor: ",corEstimate))

  }

dev.off()

integrated_4497<-integrated[["4497"]]

integrated_4497 <- ScaleData(integrated_4497, verbose = FALSE)
integrated_4497 <- RunPCA(integrated_4497, npcs = 100, verbose = FALSE)
integrated_4497 <- RunUMAP(integrated_4497, reduction = "pca", dims = 1:100)
integrated_4497 <- FindNeighbors(integrated_4497, dims.use = 1:100, resolution = 0.8)
integrated_4497 <- FindClusters(integrated_4497, dims.use = 1:100, resolution = 0.8)
integrated_4497 <- RunTSNE(integrated_4497, dims.use = 1:100,check_duplicates=FALSE)



pdf(paste0(figDir,"CITE_integrated_4497alone_tsne.pdf"),width=12,height=8)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F","darkblue","orange")
DimPlot(object = integrated_4497, 
         label = T, 
         pt.size = 2, 
         label.size = 6,
         cols = cols
)
dev.off()



genes<-c("nFeature_RNA","percent.mt",paste0(c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,"CITE_integrated_4497alone_feature_genes.pdf"),width=12,height=16)
FeaturePlot(object = integrated_4497, features = genes, ncol = 3, pt.size = 0.5,cols = c("grey", "blue"), 
            reduction = "tsne")

dev.off()
#}



integrated_4497[["ADT"]] <- CreateAssayObject(counts = predictions1@data)
integrated_4497 <- NormalizeData(integrated_4497, assay = "ADT", normalization.method = "CLR")
integrated_4497 <- ScaleData(object = integrated_4497, assay = "ADT")


pdf(paste0(figDir,"CITE_integrated_4497alone_feature_ADT.pdf"),width=12,height=16)
FeaturePlot(object = integrated_4497, features = row.names(integrated_4497$ADT@data), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)
dev.off()

cite1S<-cite1[row.names(cite1) %in% row.names(cite2),]
integrated_4497[["ADT"]] <- CreateAssayObject(counts = cite1S)
integrated_4497 <- NormalizeData(integrated_4497, assay = "ADT", normalization.method = "CLR")
integrated_4497 <- ScaleData(object = integrated_4497, assay = "ADT")


pdf(paste0(figDir,"CITE_original_4497alone_feature_ADT.pdf"),width=12,height=16)
FeaturePlot(object = integrated_4497, features = row.names(integrated_4497$ADT@data), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)
dev.off()




