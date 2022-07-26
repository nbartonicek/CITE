library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(sctransform)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/QC/")




projectnames=paste0("4T1_3x.",1:10)
sampleName="4T1_3x"
treatment="v3"

countsTotal<-list()
lossTotal<-list()
abCounts<-list()
cellNums<-list()
for(projectname in projectnames){
	cat(projectname)
	cat("\n")
	treatment="v3"
	#treatment="mt0.5_nGene7000"
	system(paste0("mkdir -p ",figDir))
	citeDir=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,"_CITE.v3")

  ############ Loss plot ###########

  #Import data from CITE-seq-count report
  #Raw:

	citeReport<-paste0(citeDir,"/run_report.yaml")
	citeReportDF<-read.table(citeReport,sep=":")
	raw_count<-as.numeric(as.character(citeReportDF[4,"V2"]))


	#Assigned to barcode white list vs unmapped:
	matrix_dir = paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,"_CITE.v3/read_count/")
	barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
	features.path <- paste0(matrix_dir, "features.tsv.gz")
	matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
	mat <- readMM(file = matrix.path)
	feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
	barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
	colnames(mat) = barcode.names$V1
	rownames(mat) = feature.names$V1
	cite<-as.data.frame(as.matrix(mat))
	rownames(cite) <- paste0("CITE_", rownames(cite))
	rownames(cite) <- gsub("-\\w{15}$","",rownames(cite))
	cite <-cite[,order(colnames(cite))]

	totalReads<-sum(cite)
	unmappedReads<-rowSums(cite)["CITE_unmapped"]
	barcodedReads<-totalReads-unmappedReads

	#UMIs
	matrix_dir = paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,"_CITE.v3/umi_count/")
	barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
	features.path <- paste0(matrix_dir, "features.tsv.gz")
	matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
	mat <- readMM(file = matrix.path)
	feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
	barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
	colnames(mat) = barcode.names$V1
	rownames(mat) = feature.names$V1
	cite<-as.data.frame(as.matrix(mat))
	rownames(cite) <- paste0("CITE_", rownames(cite))
	rownames(cite) <- gsub("-\\w{15}$","",rownames(cite))
	cite <-cite[,order(colnames(cite))]


	cellNums[[projectname]]<-dim(cite)[2]


	totalUMIs<-sum(cite)
	unmappedUMIs<-rowSums(cite)["CITE_unmapped"]
	mappedUMIs<-totalUMIs-unmappedUMIs

	top5ABs<-rowSums(cite)
	top5ABs<-top5ABs[!names(top5ABs) %in% "CITE_unmapped"]
	abCounts[[projectname]]<-length(top5ABs[top5ABs>1000])
	cat(length(top5ABs[top5ABs>1000]))
	cat("\n")
	top5ABs<-sort(top5ABs,decreasing=T)[1:5]

########### counts #############
	
  counts<-list()
  counts[["raw"]]<-raw_count
  counts[["barcoded_reads"]]<-totalReads
  #counts[["unmapped_reads"]]<-unmappedReads
  counts[["AB_reads"]]<-totalReads-unmappedReads
  counts[["AB_UMIs"]]<-mappedUMIs
  #counts[["unmapped_UMIs"]]<-unmappedUMIs
  counts[["top5_ABs_UMIs"]]<-sum(top5ABs)

  countsTotal[[projectname]]<-unlist(counts)

}

df<-do.call("rbind",countsTotal)
df<-as.data.frame(df)
df$saturation<-1-df$AB_UMIs.CITE_unmapped/df$AB_reads.CITE_unmapped

pdf(paste0(figDir,"4T1_3x.saturation.pdf"),width=6,height=6)
plot(1:10,df$saturation)
dev.off()


projectnames="4497-2"
sampleName="4497-2_CELLSUS"
treatment="v3"

finalObject<-paste0("../project_results/Robjects/",sampleName,"_",treatment,"_CITE.Rdata")
load(finalObject)
citeNorm<-citeCells@assays$ADT@data

citeNormShort<-citeNorm[c(1:12,89:96),]

data=t(citeNormShort)
dataM<-melt(data)
colnames(dataM)<-c("cell","AB","count")

pdf(paste0(figDir,"allABs_4497-2.pdf"),width=12,height=8)
p<-ggplot(dataM,aes(count,group=AB))
p<-p+geom_bar(width = 0.04)
p<-p+facet_wrap(~AB,nrow=4)
p<-p+theme_minimal()
p
dev.off()

dataM$logCount<-log1p(dataM$count)

pdf(paste0(figDir,"allABs_4497-2.log.pdf"),width=12,height=8)
p<-ggplot(dataM,aes(logCount,group=AB))
p<-p+geom_bar()
p<-p+facet_wrap(~AB,nrow=4)
p
dev.off()




















