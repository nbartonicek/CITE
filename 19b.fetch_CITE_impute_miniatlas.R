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
