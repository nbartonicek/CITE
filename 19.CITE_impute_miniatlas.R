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


#load in the total miniatlas
homeDir="../../../projects/CITE/"
RobjectsDir=paste0(homeDir,"project_results/Robjects/")
finalObject<-paste0(RobjectsDir,"03_seurat_CCA_aligned_processed.Rdata")
figDir<-paste0(homeDir,"project_results/figures/imputed_miniatlas/")
system(paste("mkdir -p",figDir))
#citeCells<-readRDS(finalObject)
#totalbcs<-colnames(citeCells)

#sampleNames<-unique(gsub("_.*","",totalbcs))
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
	outWhitelist<-paste0("../raw_files/",shortName,"/CITE_fastqs/",shortName,"_whitelist_miniatlas.txt")
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

if(!file.exists(ABsListFile)){
	for(sampleName in gsub("CID","",SSsampleNames_CITE)){
		cat(sampleName)
		cat("\n")
		load(paste0(RobjectsDir,sampleName,"_v3_CITE.Rdata"))
		#cat(a)
		#cat("\n")
		cite<-cells@assays$ADT@data
		citeRaw<-cells@assays$ADT@counts


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
#then for the normalised
workingCITEMelt<-lapply(workingCITENorm,function(x){melt(as.matrix(x))})
workingCITENarrow<-do.call(rbind,workingCITEMelt)
workingCITEWideNorm<-acast(workingCITENarrow, Var1 ~ Var2)
workingCITEWideNorm[is.na(workingCITEWideNorm)]=0

#impute
cellTypes<-c("B_cells","Myeloid","STROMAL","T_cells")
imputedSeurat4celltypes<-list()
imputedSeurat4celltypesFile<-paste0(RobjectsDir,"imputedSeurat4celltypes.Rdata")
imputedSeurat4celltypesNorm<-list()
imputedSeurat4celltypesNormFile<-paste0(RobjectsDir,"imputedSeurat4celltypesNorm.Rdata")


imputeCITE<-function(cells,cellType,norm){
	cite<-workingCITEWide[,colnames(workingCITEWide) %in% colnames(cells)]
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

imputedCellsMiniatlasFile<-paste0(RobjectsDir,"imputed_alltypes_raw.Rdata")
result<-list()
if(!file.exists(imputedCellsMiniatlasFile)){
	for(cellType in cellTypes){
		cat(cellType)
		cat("\n")
		finalObject<-paste0(RobjectsDir,"RDATA_02_PROCESSED_FILTERED_",cellType,".Rdata")
		citeCells<-readRDS(finalObject)
		for(norm in c(FALSE,TRUE)){
			cat(norm)
			cat("\n")
			cells<-imputeCITE(citeCells,cellType,norm)
			type<-ifelse(norm==TRUE,"norm","raw")
			result[[paste0(cellType,"_",type)]]<-cells
			save(cells,file=paste0(RobjectsDir,"imputed_",cellType,"_",type,".Rdata"))
		}
	}
	save(result,file=imputedCellsMiniatlasFile)
} else {
	load(imputedCellsMiniatlasFile)
}

#check clustering on one example
cellType="Myeloid_cells"
type="raw"
#cells<-result[[paste0(cellType,"_",type)]]
load(paste0(RobjectsDir,"imputed_",cellType,"_",type,".Rdata"))
reductions=c("TSNESIG10","TSNESIG20","TSNECOMBINED30","TSNESIG20")
names(reductions)<-cellTypes
reductionsUmap=c("UMAPSIG10","UMAPSIG20","UMAPCOMBINED30","UMAPSIG20")
names(reductionsUmap)<-cellTypes

#then for each cell type and normalization overlay CITE
for(cellType in cellTypes){
	cat(cellType)
	cat("\n")
	outDir<-paste0(figDir,cellType,"/")
	system(paste("mkdir -p",outDir))
	for(norm in c(FALSE,TRUE)){
		type<-ifelse(norm==TRUE,"norm","raw")
		load(paste0(RobjectsDir,"imputed_",cellType,"_",type,".Rdata"))

#find markers

		adt.markers <- FindAllMarkers(cells, assay = "ADT", only.pos = TRUE)

		deABs<-unique(adt.markers$gene)
		deABs<-sort(deABs)
		pdf(paste0(outDir,"imputed_CITE_clusters_",cellType,"_",type,".pdf"),width=12,height=20)
		p<-DoHeatmap(cells, features = unique(adt.markers$gene),raster=F, assay = "ADT", angle = 90) + NoLegend()
		print(p)
		dev.off()
		
		#plot all cite data
		nABs<-length(deABs)
		nPlots<- nABs %/% 20
		if(nABs %% 20 > 0){nPlots<-nPlots+1}
		for(i in 1:nPlots){
			range1=(i-1)*20+1
			range2=(i-1)*20+20
			pdf(paste0(outDir,"feature_plot_cite_",paste0(cellType,"_",type),"_",i,".pdf"),width=12,height=10)
			p<-FeaturePlot(cells, features = deABs[range1:range2], min.cutoff = "q05", max.cutoff = "q95", 
			    ncol = 4, cols = c("lightgrey", "blue"), reduction = reductionsUmap[cellType], pt.size = 0.5)
			print(p)
			dev.off()
		}
	}
}
		
#replot based on CITE clustering

for(cellType in cellTypes){
	cat(cellType)
	cat("\n")
	outDir<-paste0(figDir,cellType,"/")
	system(paste("mkdir -p",outDir))
	for(norm in c(FALSE,TRUE)){
		type<-ifelse(norm==TRUE,"norm","raw")
		load(paste0(RobjectsDir,"imputed_",cellType,"_",type,".Rdata"))

		adt.markers <- FindAllMarkers(cells, assay = "ADT", only.pos = TRUE)
		deABs<-unique(adt.markers$gene)
		deABs<-sort(deABs)


		DefaultAssay(cells) <- "ADT"
		#cells <- RunPCA(cells, features = row.names(cells@assays$ADT@data), reduction.name = "pca_adt", reduction.key = "pca_adt_", 
		#    verbose = FALSE)
		#)
		cells <- RunTSNE(cells, features = deABs, reduction.name = "tsne_adt", reduction.key = "tsne_adt_", 
		    verbose = FALSE)
		#plot all cite data
		nABs<-length(deABs)
		nPlots<- nABs %/% 20
		if(nABs %% 20 > 0){nPlots<-nPlots+1}
		for(i in 1:nPlots){
			range1=(i-1)*20+1
			range2=(i-1)*20+20
			pdf(paste0(outDir,"feature_plot_CITE-based_",paste0(cellType,"_",type),"_",i,".pdf"),width=12,height=10)
			p<-FeaturePlot(cells, features = deABs[range1:range2], min.cutoff = "q05", max.cutoff = "q95", 
			    ncol = 4, cols = c("lightgrey", "blue"), reduction = "tsne_adt", pt.size = 0.5)
			print(p)
			dev.off()
		}
	}
}













