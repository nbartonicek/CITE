library(Seurat)
library(dplyr)
library(Matrix)
#library(devtools)
library(sctransform)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(plyr)

library(R.utils)
#for testing

args <- R.utils::commandArgs(asValues=TRUE)
if (!is.null(args[["eliminateSample"]])){eliminateSample = args$eliminateSample} 

#load in the total miniatlas
homeDir="/share/ScratchGeneral/nenbar/projects/CITE/"
RobjectsDir=paste0(homeDir,"project_results/Robjects/")
finalObject<-paste0(RobjectsDir,"03_seurat_CCA_aligned_processed.Rdata")
figDir<-paste0(homeDir,"project_results/figures/imputed_miniatlas/")
system(paste("mkdir -p",figDir))
#citeCells<-readRDS(finalObject)
#totalbcs<-colnames(citeCells)


#load in the subsections
cellTypes<-c("B_cells","Myeloid","STROMAL","T_cells")
ssBCs<-list()
ssBCsFile<-paste0(RobjectsDir,"ssBCsFile.Rdata")
if(!file.exists(ssBCsFile)){
	for(cellType in cellTypes){
		cat(cellType)
		cat("\n")
		finalObject<-paste0(RobjectsDir,"RDATA_02_PROCESSED_FILTERED_",cellType,".Rdata")
		citeCells<-readRDS(finalObject)
		ssBCs[[cellType]]<-colnames(citeCells)
	}
	save(ssBCs,file=ssBCsFile)
} else {load(ssBCsFile)}

subsampleBCs<-unlist(ssBCs)
SSsampleNames<-unique(gsub("_.*","",subsampleBCs))

#select sampleNames with CITE
SSsampleNames_CITE<-SSsampleNames[grepl("3838|3946|4040|4378|4515",SSsampleNames)]

#create individual whitelists from SSsampleName
for(sampleName in SSsampleNames_CITE){
	cat(sampleName)
	cat("\n")
	#get all barcodes
	bcs<-subsampleBCs[grepl(sampleName,subsampleBCs)]
	bcs<-gsub(".*_","",bcs)
	cat(length(bcs))
	cat("\n")
	
	#export them
	shortName<-gsub("CID","",sampleName)
	outWhitelist<-paste0(homeDir,"raw_files/",shortName,"/CITE_fastqs/",shortName,"_whitelist_miniatlas.txt")
	write(bcs,outWhitelist)
}

#for each sample calculate new counts for subsetted data only
#How to select which antibody will be imputed if we do not know signal/noise
#Take only those that are usable markers (DE in clusters) to select biological clusters
#Or take those 
#and then selecting those that were biomarkers for the cell type


ABsList<-list()
ABsListNorm<-list()
ABsListFile<-paste0(RobjectsDir,"ABsList.Rdata")
ABsListNormFile<-paste0(RobjectsDir,"ABsListNorm.Rdata")

#for each sample extract data
workingCITE<-list()
workingCITEFile<-paste0(RobjectsDir,"workingCITE.Rdata")
workingCITENorm<-list()
workingCITENormFile<-paste0(RobjectsDir,"workingCITENorm.Rdata")

rawCiteCounts<-list()
if(!file.exists(ABsListFile)){
	for(sampleName in gsub("CID","",SSsampleNames_CITE)){
		cat(sampleName)
		cat("\n")
		load(paste0(RobjectsDir,sampleName,"_v3_CITE.Rdata"))
		#cat(a)
		#cat("\n")
		cite<-cells@assays$ADT@data
		citeRaw<-cells@assays$ADT@counts
		rawCiteCounts[[sampleName]]<-citeRaw[rowSums(citeRaw)>1000,]


		#decide which antibody will be used by loading the Robjects and selecting 

		adt.markers <- FindAllMarkers(cells, assay = "ADT", only.pos = TRUE,verbose=F)
		adt.markers <- adt.markers[!grepl("unmapped",adt.markers$gene),]
	 	cat("Total marker antibodies:\t")
		cat(length(unique(adt.markers$gene)))
		cat("\n")
		ABsList[[sampleName]]<-unique(adt.markers$gene)

		igGs<-cite[grepl("IgG\\d",row.names(cite)),]
		igGsShort<-igGs[,colSums(igGs)>0]
		meanIgG<-apply(igGsShort,2,mean)
		citeShort<-cite[,colSums(igGs)>0]
		citeShort<-as.data.frame(citeShort)

		citeAllNorm<-t(apply(citeShort,1,function(x){as.integer(x/meanIgG)}))
		colnames(citeAllNorm)<-colnames(citeShort)
		citeAllNorm<-as.data.frame(citeAllNorm)


		citeAllNorm$id<-row.names(citeAllNorm)
		allBCs<-colnames(cells@assays$RNA@data)

		citeAll<-as.data.frame(matrix(0,dim(citeAllNorm)[1],length(allBCs)))
		row.names(citeAll)<-row.names(citeAllNorm)
		colnames(citeAll)<-allBCs

		citeAll$id<-row.names(citeAll)
		citeL<-list(cite=citeAllNorm,citeAll=citeAll)

		matrix.df <- ldply(citeL, melt)
		sum.matrix <- acast(matrix.df, id ~ variable, sum)

		citeAll<-sum.matrix



		cells[["ADT"]] <- CreateAssayObject(counts = citeAll)
		cells <- NormalizeData(cells, assay = "ADT", normalization.method = "CLR")
		cells <- ScaleData(object = cells, assay = "ADT")


		#for each 
		adt.markers <- FindAllMarkers(cells, assay = "ADT", only.pos = TRUE,verbose=F)
		adt.markers <- adt.markers[!grepl("unmapped",adt.markers$gene),]
	 	cat("Total marker normalized antibodies:\t")
		cat(length(unique(adt.markers$gene)))
		cat("\n")
		ABsListNorm[[sampleName]]<-unique(adt.markers$gene)

		#take antibodies that worked
		liveABs<-ABsList[[sampleName]]
		#eliminate IgG
		liveABs<-liveABs[!grepl("IgG",liveABs)]
		#add sample IDs
		#cite<-cells@assays$ADT@counts
		colnames(citeRaw)<-paste0("CID",sampleName,"_",colnames(citeRaw))
		workingCITE[[sampleName]]<-citeRaw[row.names(citeRaw) %in% liveABs,]

		#repeat for IgG normalised
		liveABs<-ABsListNorm[[sampleName]]
		#eliminate IgG
		liveABs<-liveABs[!grepl("IgG",liveABs)]
		colnames(citeAll)<-paste0("CID",sampleName,"_",colnames(citeAll))
		
		workingCITENorm[[sampleName]]<-citeAll[row.names(citeAll) %in% liveABs,]

	}
	save(ABsList,file=ABsListFile)
	save(ABsListNorm,file=ABsListNormFile)
	save(workingCITE,file=workingCITEFile)
	save(workingCITENorm,file=workingCITENormFile)

} else {
	load(ABsListFile);
	load(ABsListNormFile);
	load(workingCITEFile);
	load(workingCITENormFile)
}

workingCITEMelt<-lapply(workingCITE,function(x){melt(as.matrix(x))})
workingCITENarrow<-do.call(rbind,workingCITEMelt)
workingCITEWide<-acast(workingCITENarrow, Var1 ~ Var2)
workingCITEWide[is.na(workingCITEWide)]=0

#impute
#cellTypes<-c("T_cells","B_cells","Myeloid","STROMAL",)
cellTypes<-c("T_cells")
imputedSeurat4celltypes<-list()
imputedSeurat4celltypesFile<-paste0(RobjectsDir,"imputedSeurat4celltypes.Rdata")


imputeCITE<-function(cells,cellType,norm,eliminate){
	cite<-workingCITEWide[,colnames(workingCITEWide) %in% colnames(cells)]
	#do not impute with all samples, eliminate one in turn
	cite<-cite[,!grepl(eliminate,colnames(cite))]
	
	if(norm){
		cite<-workingCITEWideNorm[,colnames(workingCITEWideNorm) %in% colnames(cells)]
	}
	cite<-as.data.frame(cite)
	#fix missing zeroes
	cite$id<-row.names(cite)
	allBCs<-colnames(cells@assays$RNA@data)

	citeAll<-as.data.frame(matrix(0,dim(cite)[1],length(allBCs)))
	row.names(citeAll)<-row.names(cite)
	colnames(citeAll)<-allBCs

	citeAll$id<-row.names(citeAll)
	citeL<-list(cite=cite,citeAll=citeAll)

	matrix.df <- ldply(citeL, melt)
	#some values are NA, so turn them into 0
	matrix.df$value[is.na(matrix.df$value)]<-0
	sum.matrix <- acast(matrix.df, id ~ variable, sum)

	citeAll<-sum.matrix
	citeAll<-as.data.frame(citeAll)
	citeAll<-citeAll[colnames(citeAll) %in% colnames(cells),]
	cells[["ADT"]] <- CreateAssayObject(counts = citeAll)
	cells <- NormalizeData(cells, assay = "ADT", normalization.method = "CLR")
	cells <- ScaleData(object = cells, assay = "ADT")


	allCells<-colnames(cells@assays$RNA@data)
	nonQuery<-allCells[!(allCells %in% colnames(cite))]
	reference <- subset(cells,cells=colnames(cite) )
	query <- subset(cells,cells=nonQuery )

	query.anchors <- FindTransferAnchors(reference = reference, query = query, project.query=T)
	predictions1 <- TransferData(anchorset = query.anchors, refdata = reference@assays$ADT@counts)
	predictionsDF<-as.data.frame(predictions1@data)
	cite<-cite[order(row.names(cite)),]
	predictionsDF<-predictionsDF[row.names(predictionsDF) %in% row.names(cite),]
	predictionsDF<-predictionsDF[order(row.names(predictionsDF)),]

	cite<-cbind(cite,predictionsDF)
	cite<-cite[,colnames(cite) %in% colnames(cells@assays$RNA@data)]


	cells[["ADT"]] <- CreateAssayObject(counts = cite)
	cells <- NormalizeData(cells, assay = "ADT", normalization.method = "CLR")
	cells <- ScaleData(object = cells, assay = "ADT")

	
}



#for(SSsampleName_CITE in eliminateSample){
#	imputedCellsMiniatlasFile<-paste0(RobjectsDir,"imputed_alltypes_without_",SSsampleName_CITE,".Rdata")
#	result<-list()
#	if(!file.exists(imputedCellsMiniatlasFile)){
#		for(cellType in cellTypes){
#			cat(cellType)
#			cat("\n")
#			finalObject<-paste0(RobjectsDir,"RDATA_02_PROCESSED_FILTERED_",cellType,".Rdata")
#			citeCells<-readRDS(finalObject)
#			for(norm in c(FALSE)){
#				cat(norm)
#				cat("\n")
#				cells<-imputeCITE(citeCells,cellType,norm,SSsampleName_CITE)
#				type<-ifelse(norm==TRUE,"norm","raw")
#				result[[paste0(cellType,"_",type)]]<-cells
#				save(cells,file=paste0(RobjectsDir,"imputed_",cellType,"_without_",SSsampleName_CITE,".Rdata"))
#			}
#		}
#		save(result,file=imputedCellsMiniatlasFile)
#	} else {
#		load(imputedCellsMiniatlasFile)
#	}
#}



#check clustering on one example
#cellType="STROMAL"
#type="raw"
#load(paste0(RobjectsDir,"imputed_",cellType,"_",type,".Rdata"))
reductions=c("TSNESIG10","TSNESIG20","TSNECOMBINED30","TSNESIG20")
names(reductions)<-cellTypes
reductionsUmap=c("UMAPSIG10","UMAPSIG20","UMAPCOMBINED30","UMAPSIG20")
names(reductionsUmap)<-cellTypes

cellType="T_cells"
SSsampleName_CITE="CID4515"

#1. load imputed object without 4515
inFile=paste0(RobjectsDir,"imputed_",cellType,"_without_",SSsampleName_CITE,".Rdata")
load(inFile)
cellsImputedWithoutSample=cells
#2. load imputed object with 4515
type="raw"
inFileTotal=paste0(RobjectsDir,"miniatlas_imputed/imputed_",cellType,"_",type,".Rdata")
load(inFileTotal)

#3. plot for on 4515 cells only
#get 4515 cells
cellIDs<-colnames(cellsImputedWithoutSample)
cellIDs<-cellIDs[grepl(SSsampleName_CITE,cellIDs)]
cellsImputedWithoutSample_ss<-subset(cellsImputedWithoutSample,cells=cellIDs)
cells_ss<-subset(cells,cells=cellIDs)
sum(colnames(cellsImputedWithoutSample_ss) %in% colnames(cells_ss))

#4. make the umap for the cells in mind and plot featureplots and CITE over them
cite4515file<-paste0(RobjectsDir,"only4515_miniatlas_umap.Rdata")

if(!file.exists(cite4515file)){
	npcs=100
	citeCells<- cells_ss
	citeCells <- SCTransform(object = citeCells, verbose = FALSE)
	citeCells <- RunPCA(object = citeCells, verbose = FALSE,npcs = npcs)
	citeCells <- JackStraw(object = citeCells, dims=100,num.replicate=100,verbose=T)
	citeCells <- ScoreJackStraw(citeCells,dims=1:100)
	pcsData<-as.data.frame(citeCells@reductions$pca@jackstraw@overall.p.values)
	pcs<-pcsData[pcsData$Score<0.05,"PC"]
	citeCells <- FindNeighbors(object = citeCells, dims = pcs, verbose = FALSE)
	citeCells <- FindClusters(citeCells, dims.use = pcs, resolution = 0.6)
	citeCells <- RunUMAP(citeCells, dims = pcs)
	save(citeCells,file=cite4515file)
} else {load(cite4515file)}

genes<-c("nFeature_RNA","percent.mt",paste0(c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,cellType,"_",SSsampleName_CITE,"_feature_vst.pdf"),width=12,height=16)
FeaturePlot(object = citeCells, features = genes, ncol = 3, pt.size = 0.5,cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()

abs<-row.names(citeCells@assays$ADT@data)

selected_abs<-c("CITE-CD4","CITE-CD8a","CITE-CD16","CITE-CD57","CITE-PD-1","CITE-PD-L1","CITE-PD-L2","CITE-CD103")
pdf(paste0(figDir,"feature_plot_cite_",cellType,"_",SSsampleName_CITE,"_selected_ABs.pdf"),width=12,height=4)
FeaturePlot(citeCells, features = selected_abs, min.cutoff = "q05", max.cutoff = 2, 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()

pdf(paste0(figDir,"feature_plot_cite_",cellType,"_",SSsampleName_CITE,"_1.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = abs[1:20], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()

pdf(paste0(figDir,"feature_plot_cite_",cellType,"_",SSsampleName_CITE,"_1.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = abs[1:20], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",cellType,"_",SSsampleName_CITE,"_2.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = abs[21:40], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",cellType,"_",SSsampleName_CITE,"_3.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = abs[41:60], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",cellType,"_",SSsampleName_CITE,"_4.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = abs[61:80], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",cellType,"_",SSsampleName_CITE,"_5.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = abs[81:100], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",cellType,"_",SSsampleName_CITE,"_6.pdf"),width=12,height=5)
FeaturePlot(citeCells, features = abs[101:length(abs)[1]], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()



#3a. plot for on 4515 cells only
#get 4515 cells





cells_ws<-subset(cellsImputedWithoutSample,cells=cellIDs)

adt.markers <- FindAllMarkers(cells_ws, assay = "ADT", only.pos = TRUE)
adtL<-split(adt.markers,adt.markers$gene)
minVal<-lapply(adtL,function(x){min(x[,1])})
minVal<-unlist(minVal)

deABs<-names(minVal)[order(minVal)]
imputedCITE<-cells_ws@assays$ADT@counts
imputedCITE<-imputedCITE[row.names(imputedCITE) %in% deABs,]
imputedCITE<-imputedCITE[deABs,]
row.names(imputedCITE)<-paste0(row.names(imputedCITE),"_imputed")

cite<-citeCells@assays$ADT@counts
cite<-cite[row.names(cite) %in% deABs,]
cite<-cite[deABs,]
cite<-rbind(cite,imputedCITE)

abs<-rep(deABs,each=2)
abs<-paste0(abs,c("","_imputed"))
cite<-cite[abs,]
abs<-gsub("_","-",abs)
citeCells[["ADT"]] <- CreateAssayObject(counts = cite)
citeCells <- NormalizeData(citeCells, assay = "ADT", normalization.method = "CLR")
citeCells <- ScaleData(object = citeCells, assay = "ADT")



pdf(paste0(figDir,"imputed_significant_",cellType,"_",SSsampleName_CITE,"_1.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = abs[1:16], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"imputed_significant_",cellType,"_",SSsampleName_CITE,"_2.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = abs[17:32], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"imputed_significant_",cellType,"_",SSsampleName_CITE,"_3.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = abs[33:48], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"imputed_significant_",cellType,"_",SSsampleName_CITE,"_4.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = abs[49:64], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"imputed_significant_",cellType,"_",SSsampleName_CITE,"_5.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = abs[65:80], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"imputed_significant_",cellType,"_",SSsampleName_CITE,"_6.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = abs[81:96], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()

pdf(paste0(figDir,"imputed_significant_",cellType,"_",SSsampleName_CITE,"_7.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = abs[97:112], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
dev.off()

################################
# 1. first add correlation coefficients

#for each antibody that exists in both (?) add a correlation coefficient
res<-list()
citeNorm<-citeCells@assays$ADT@scale.data
for(i in 1:length(deABs)){
	res[[deABs[i]]]<-cor.test(as.numeric(citeNorm[2*i-1,]),as.numeric(citeNorm[2*i,]))$estimate
}

#corre

# then add correlation coefficients of each with RNA
# one plot 
#CITE-CD19.cor 
#       0.0595222833        0.0972415956        0.1114999610        0.2145816794 
#      CITE-CD86.cor     CITE-CD45RA.cor      CITE-TIGIT.cor        CITE-CD3.cor 
#       0.2435019699        0.2496093117        0.3009984685        0.4441411345 
#      CITE-CD8a.cor       CITE-CD16.cor        CITE-CD4.cor       CITE-PD-1.cor 
#       0.4753007218        0.5795148485        0.6594383225        0.6678998062






