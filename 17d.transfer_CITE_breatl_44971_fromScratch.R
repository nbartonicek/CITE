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

homeDir="/share/ScratchGeneral/nenbar/projects/CITE/"
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
finalObject<-paste0(RobjectsDir,"03_seurat_CCA_aligned_processed.Rdata")
#finalObject<-"/paella/TumourProgressionGroupTemp/projects/brca_mini_atlas/analysis/randomly_subsetted_20k_cells_for_testing/sDataRandomSubset.Rdata"
citeCells<-readRDS(finalObject)
citeCells1<-citeCells



######## integrate ########
citeCells1Counts<-citeCells1@assays$RNA@counts
citeCells2Counts<-citeCellsShort@assays$RNA@counts
citeCells1Counts<-citeCells1Counts[row.names(citeCells1Counts) %in% row.names(citeCells2Counts),]
citeCells2Counts<-citeCells2Counts[row.names(citeCells2Counts) %in% row.names(citeCells1Counts),]

c1<-CreateSeuratObject(citeCells1Counts)
c2<-CreateSeuratObject(citeCells2Counts)

integratedFile=paste0("../project_results/Robjects/CITE_integrated_transfer_4497_miniatlas.Rdata")
if(!file.exists(integratedFile)){
	integrated<-list()

	integrated[["integrated"]]<-NormalizeData(c1,verbose=FALSE)
	integrated[["4497"]]<-NormalizeData(c2,verbose=FALSE)
	integrated[["integrated"]] <- FindVariableFeatures(integrated[["integrated"]], selection.method = "vst", nfeatures = 2000, 
	        verbose = FALSE)
	integrated[["4497"]] <- FindVariableFeatures(integrated[["4497"]], selection.method = "vst", nfeatures = 2000, 
	        verbose = FALSE)
	integrated[["4497"]][["ADT"]] <- CreateAssayObject(counts = cite)
	integrated[["4497"]] <- NormalizeData(integrated[["4497"]], assay = "ADT", normalization.method = "CLR")
	integrated[["4497"]] <- ScaleData(object = integrated[["4497"]], assay = "ADT")
} else {load(integratedFile)}

anchorsFile<-paste0("../project_results/Robjects/CITE_integrated_transfer_4497_miniatlas.postIntegration.Rdata")

if(!file.exists(anchorsFile)){

	anchors <- FindIntegrationAnchors(object.list = integrated, dims = 1:30)
	cite.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
} else {load(anchorsFile)}


DefaultAssay(cite.integrated) <- "integrated"
cite.integrated <- ScaleData(cite.integrated, verbose = FALSE)
cite.integrated <- RunPCA(cite.integrated, npcs = 30, verbose = FALSE)

cite.integrated@meta.data$orig.ident=rep(c("integrated","4497"),times=c(dim(citeCells1Counts)[2],dim(citeCells2Counts)[2]))

query.anchorsFile=paste0("../project_results/Robjects/CITE_integrated_transfer_4497_miniatlas.query.anchors3.Rdata")
if(!file.exists(query.anchorsFile)){
	query <- integrated[["integrated"]]
	query.anchors <- FindTransferAnchors(reference = cite.integrated, query = query, 
    dims = 1:30, project.query=T)
	save(query.anchors,file=query.anchorsFile)
} else {load(query.anchorsFile)}

predictions1 <- TransferData(anchorset = query.anchors, refdata = cite.integrated@assays$ADT@counts, 
    dims = 1:30)
citeNew<-as.data.frame(predictions1@data)
save(citeNew,file=paste0("../project_results/Robjects/44971-imputed_miniatlas_CITE.Rdata"))


pdf(paste0(figDir,"CITE_imputed_miniatlas_feature_ADT_all.pdf"),width=12,height=16)
FeaturePlot(object = citeCells1, features = row.names(citeCells1$ADT@data), reduction="UMAPC",min.cutoff = "q05", max.cutoff = "q95", ncol = 4)
dev.off()





