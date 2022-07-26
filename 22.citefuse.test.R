library(CiteFuse)
library(scater)
library(SingleCellExperiment)
library(DT)
library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)

#library(DropletUtils)
library(R.utils)
library(plyr)
library(reshape2)
library(ggpubr)


#questions to answer for tomorrow
#1. how does citefuse work
#2. what is spectral clustering
#3. what is importance score
#4. can it cluster on both

#1. results across tissues
#a. counts
#b. normalised counts
#c. importance scores



#what is spectral clustering

#is the clustering done on one or both?

###run citefuse on all samples 

###save object

###compare the number of usable "significant" antibodies


#for testing
fileName="RDATA_01_PROCESSED_FILTERED_IMPUTEDCITE_RAW_tcells.Rdata"

args <- R.utils::commandArgs(asValues=TRUE)
if (!is.null(args[["fileName"]])){fileName = args$fileName} 

projectname="CITE"
homeDir=paste0("/share/ScratchGeneral/nenbar/projects/",projectname,"/")

treatment="v4"
#Treatment V4 includes the following cutoffs:
#	min.cells = 3, 
#	min.features = 100
#	nFeature_RNA > 250 
#	percent.mt < 25
#	npcs=100
min.cells = 3
min.features = 100
nFeature_RNA = 250 
percent.mt = 25
npcs=100

sampleName=gsub(paste0(".Rdata"),"",fileName)
sampleName=gsub(paste0(".*_RAW"),"",sampleName)
sampleName=gsub(paste0(".*_"),"",sampleName)

figDir=paste0(homeDir,"project_results/figures/citefuse/")
system(paste0("mkdir -p ",figDir))

RobjectsDir=paste0(homeDir,"project_results/Robjects/citefuse/")
system(paste0("mkdir -p ",RobjectsDir))
tableDir=paste0(homeDir,"project_results/tables/citefuse/")
system(paste0("mkdir -p ",tableDir))


#load data from R object
hashedFile=paste0(RobjectsDir,fileName)
citeCells<-readRDS(hashedFile)
#samplename="log"
#make a citefuse object if it doesn't exist
fuseFile=paste0(RobjectsDir,sampleName,".Rdata")

if(!file.exists(fuseFile)){

	cells<-list()
	cells[["RNA"]]<-citeCells@assays$RNA@counts
	#cite<-as.matrix(citeCells@assays$ADT@data)
	#citeL<-apply(cite,1,function(x){100^x-1})
	#citeDF<-t(citeL)

	#cells[["ADT"]]<-citeDF
	cells[["ADT"]]<-citeCells@assays$ADT@data


	sce_citeseq <- preprocessing(cells)
	sce_citeseq <- scater::logNormCounts(sce_citeseq)
	#sce_citeseq <- normaliseExprs(sce_citeseq, altExp_name = "ADT", transform = "log")
	system.time(sce_citeseq <- CiteFuse(sce_citeseq))

	SNF_W_clust <- spectralClustering(metadata(sce_citeseq)[["SNF_W"]], K = 20)

	pdf(paste0(figDir,sampleName,"_SNF_W_clust.pdf"),width=12,height=10)
	plot(SNF_W_clust$eigen_values)
	dev.off()

	maxVal<-as.numeric(which.max(abs(diff(SNF_W_clust$eigen_values))))
	#maxVal = 6
	SNF_W_clust <- spectralClustering(metadata(sce_citeseq)[["SNF_W"]], K = maxVal)

	sce_citeseq$SNF_W_clust <- as.factor(SNF_W_clust$labels)
	SNF_W1_clust <- spectralClustering(metadata(sce_citeseq)[["ADT_W"]], K = maxVal)

	sce_citeseq$ADT_clust <- as.factor(SNF_W1_clust$labels)
	SNF_W2_clust <- spectralClustering(metadata(sce_citeseq)[["RNA_W"]], K = maxVal)

	sce_citeseq$RNA_clust <- as.factor(SNF_W2_clust$labels)


	sce_citeseq <- reducedDimSNF(sce_citeseq,
	                             method = "UMAP", 
	                             dimNames = "UMAP_joint")
	SNF_W_louvain <- igraphClustering(sce_citeseq, method = "louvain")
	sce_citeseq$SNF_W_louvain <- as.factor(SNF_W_louvain)
	sce_citeseq <- importanceADT(sce_citeseq, 
                             group = sce_citeseq$SNF_W_louvain,
                             subsample = TRUE)

  save(sce_citeseq,file=fuseFile)
} else {load(fuseFile)}

#save importance data
imp<-metadata(sce_citeseq)[["importanceADT"]]
df<-data.frame(adt=names(imp),pval=as.numeric(imp))

write.table(df,file=paste0(tableDir,sampleName,".txt"),quote=F,sep="\t",row.names=F)

g1 <- visualiseDim(sce_citeseq, dimNames = "UMAP_joint", colour_by = "SNF_W_clust") +
  labs(title = "UMAP (SNF clustering)")
g2 <- visualiseDim(sce_citeseq, dimNames = "UMAP_joint",  colour_by = "ADT_clust") +
  labs(title = "UMAP (ADT clustering)")
g3 <- visualiseDim(sce_citeseq, dimNames = "UMAP_joint",  colour_by = "RNA_clust") +
  labs(title = "UMAP (RNA clustering)")

library(gridExtra)
pdf(paste0(figDir,sampleName,"_UMAP_clustering.pdf"),width=12,height=10)
p<-grid.arrange(g3, g2, g1, ncol = 2)
print(p)
dev.off()




pdf(paste0(figDir,sampleName,"_louvain_clustering.pdf"),width=12,height=10)
visualiseDim(sce_citeseq, dimNames = "UMAP_joint", colour_by = "SNF_W_louvain") +
  labs(title = "UMAP (SNF louvain clustering)")
dev.off()

pdf(paste0(figDir,sampleName,"_CITE_expression.pdf"),width=12,height=10)
visualiseExprs(sce_citeseq, 
               altExp_name = "ADT", 
               n = 30) +
coord_flip()
dev.off()

pdf(paste0(figDir,sampleName,"_CD4_CD8.pdf"),width=12,height=10)
visualiseExprs(sce_citeseq, altExp_name = "ADT", 
               plot = "pairwise", 
               feature_subset = c("CITE-CD4", "CITE-CD8A"))
dev.off()


pdf(paste0(figDir,sampleName,"_importance.pdf"),width=10,height=12)
visImportance(sce_citeseq, plot = "boxplot")
dev.off()

pdf(paste0(figDir,sampleName,"_importance_heatmap.pdf"),width=18,height=16)
visImportance(sce_citeseq, plot = "heatmap")
dev.off()

