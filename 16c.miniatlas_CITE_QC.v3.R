library(Seurat)
library(dplyr)
library(Matrix)
#library(devtools)
library(sctransform)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/QC_miniatlas/")
RobjectsDir=paste0(homeDir,"project_results/Robjects/")
treatment="v3"




projectnames=c("4040","3946","3838")
countsTotal<-list()
lossTotal<-list()
abCounts<-list()
cellNums<-list()
for(projectname in projectnames){
  cat(projectname)
  cat("\n")
  sampleName=paste0(projectname,"_CITE")
  treatment="v3"
  #treatment="mt0.5_nGene7000"
  system(paste0("mkdir -p ",figDir))
  citeDir=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,".v3")

  ############ Loss plot ###########

  #Import data from CITE-seq-count report
  #Raw:

  citeReport<-paste0(citeDir,"/run_report.yaml")
  citeReportDF<-read.table(citeReport,sep=":")
  raw_count<-as.numeric(as.character(citeReportDF[4,"V2"]))


  #/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_*-1*092019/run01_sunwupipeline/output/count/CID$sampleName/count_CID"$sampleName"_GRCh38_mm10/outs/raw_feature_bc_matrix/
  if(projectname=="3838"){
    inDir=paste0("/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_CITE-16092019/run01_sunwupipeline/output/count/CID",projectname,"/count_CID",projectname,"_GRCh38_mm10/outs/filtered_feature_bc_matrix/")
  } else if (projectname=="3946"){
    inDir=paste0("/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_CITE-19092019/run01_sunwupipeline/output/count/CID",projectname,"/count_CID",projectname,"_GRCh38_mm10/outs/filtered_feature_bc_matrix/")
  } else if (projectname=="4040"){
    inDir=paste0("/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_HTO-17092019/run01_sunwupipeline/output/count/CID",projectname,"/count_CID",projectname,"_GRCh38_mm10/outs/filtered_feature_bc_matrix/")
  }
  if(projectname=="3838"){
    inFile=paste0("/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_CITE-16092019/run01_sunwupipeline/output/seurat/seurat_CID",projectname,"/Output/Rdata/03_seurat_object_processed.RData")
  } else if (projectname=="3946"){
    inFile=paste0("/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_CITE-19092019/run01_sunwupipeline/output/seurat/seurat_CID",projectname,"/Output/Rdata/03_seurat_object_processed.RData")
  } else if (projectname=="4040"){
    inFile=paste0("/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_HTO-17092019/run01_sunwupipeline/output/seurat/seurat_CID",projectname,"/Output/Rdata/03_seurat_object_processed.RData")
  }
  cells<-readRDS(inFile)

  #Assigned to barcode white list vs unmapped:
  matrix_dir = paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,".v3/read_count/")
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  cite<-as.data.frame(mat)
  cite <-cite[,colnames(cite) %in% gsub(".*_","",colnames(cells$RNA@data)) ]
  rownames(cite) <- paste0("CITE_", rownames(cite))
  rownames(cite) <- gsub("-\\w{15}$","",rownames(cite))
  cite <-cite[,order(colnames(cite))]

  totalReads<-sum(cite)
  unmappedReads<-rowSums(cite)["CITE_unmapped"]
  barcodedReads<-totalReads-unmappedReads

  #UMIs
  matrix_dir = paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,".v3/umi_count/")
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  cite<-as.data.frame(mat)
  cite <-cite[,colnames(cite) %in% gsub(".*_","",colnames(cells$RNA@data)) ]
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

  df<-do.call("rbind",counts)
  colnames(df)<-"counts"
  dataM<-melt(df)
  colnames(dataM)<-c("type","class","counts")

  pdf(paste0(figDir,sampleName,"_",treatment,"_counts.pdf"),width=6,height=6)
  p<-ggplot(dataM,aes(type,counts))
  p<-p+geom_bar(stat="identity")
  p<-p+theme(axis.text.x=element_text(angle = 45, vjust = 0.5))
  print(p)
  dev.off()


  ########## loss plot ###########

  #Loss:
  #unassigned to barcodes
  #unmapped to antibodies
  #duplicated
  #absorbed_by_top5

  loss<-list()
  loss[["non_whitelist_barcode"]]<-raw_count-totalReads
  loss[["unmapped_to_AB"]]<-unmappedReads
  loss[["mapped_redundant"]]<-totalReads-unmappedReads-mappedUMIs
  loss[["absorbed_by_top5"]]<-sum(top5ABs)
  loss[["available"]]<-mappedUMIs-sum(top5ABs)

  lossTotal[[projectname]]<-unlist(loss)

  df<-do.call("rbind",loss)
  colnames(df)<-"counts"
  dataM<-melt(df)
  colnames(dataM)<-c("type","class","counts")

  pdf(paste0(figDir,sampleName,"_",treatment,"_loss.pdf"),width=6,height=6)
  p<-ggplot(dataM,aes(type,counts))
  p<-p+geom_bar(stat="identity")
  p<-p+ylim(c(0,raw_count))
  p<-p+theme(axis.text.x=element_text(angle = 45, vjust = 0.5))
  print(p)

  dev.off()

}



cols<-c("orange","red","lightblue","darkblue","purple")
cols<-rev(brewer.pal(5,"Set1"))
df<-do.call("rbind",lossTotal)
dataM<-melt(df)
colnames(dataM)<-c("sample","class","counts")
dataM$class<-gsub("\\..*","",dataM$class)
dataM$class<-factor(dataM$class,levels=names(loss))


pdf(paste0(figDir,"lossAll_counts.pdf"),width=6,height=6)
p<-ggplot(dataM,aes(sample,counts,fill=class))
p<-p+geom_bar(stat="identity")
p<-p+theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
p<-p+scale_fill_manual(values=cols)
p<-p+ylim(c(0,raw_count))
p
dev.off()


finalObject<-paste0(RobjectsDir,projectname,"_",treatment,"_CITE.Rdata")
load(finalObject)
citeNorm<-citeCells@assays$ADT@data

nCites<-dim(citeNorm)[1]
citeNorm<-citeNorm[order(rowSums(citeCells@assays$ADT@counts)),]
citeNormShort<-citeNorm[c(1:12,(nCites-7):nCites),]

data=t(citeNormShort)
dataM<-melt(data)
colnames(dataM)<-c("cell","AB","count")
dataM$AB<-factor(dataM$AB,levels=row.names(citeNormShort))
pdf(paste0(figDir,"allABs_",sampleName,".pdf"),width=10,height=10)
p<-ggplot(dataM,aes(count,group=AB))
p<-p+geom_bar()
p<-p+facet_wrap(~AB,ncol=4)
p
dev.off()


citeCounts=sort(rowSums(citeCells@assays$ADT@counts))
df=data.frame(count=citeCounts,marker=gsub("CITE-","",names(citeCounts)))

df$marker<-factor(df$marker,levels=gsub("CITE-","",names(citeCounts)))
nMarkers<-dim(df)[1]
df$group<-cut(1:nMarkers,6)

pdf(paste0(figDir,sampleName,"_reads.pdf"),width=12,height=8)
p<-ggplot(df,aes(marker,count))+geom_point()
p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
p<-p+coord_flip()
p<-p+facet_wrap(~group,,scales='free',ncol=3,drop=T)
p
dev.off()

##### saturation ######
break()


subsamples=paste0(projectname,".",1:10)
treatment="v3"

countsTotal<-list()
lossTotal<-list()
abCounts<-list()
cellNums<-list()
for(subsample in subsamples){
  cat(subsample)
  cat("\n")
  treatment="v3"
  #treatment="mt0.5_nGene7000"
  system(paste0("mkdir -p ",figDir))
  citeDir=paste0(homeDir,"raw_files/",subsample,"/CITE_fastqs/",projectname,"_CITE.v3")

  ############ Loss plot ###########

  #Import data from CITE-seq-count report
  #Raw:

  citeReport<-paste0(citeDir,"/run_report.yaml")
  citeReportDF<-read.table(citeReport,sep=":")
  raw_count<-as.numeric(as.character(citeReportDF[4,"V2"]))


  #Assigned to barcode white list vs unmapped:
  matrix_dir = paste0(citeDir,"/read_count/")
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
  matrix_dir = paste0(citeDir,"/umi_count/")
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
  cellNums[[subsample]]<-dim(cite)[2]
  totalUMIs<-sum(cite)
  unmappedUMIs<-rowSums(cite)["CITE_unmapped"]
  mappedUMIs<-totalUMIs-unmappedUMIs

########### counts #############
  
  counts<-list()
  counts[["AB_reads"]]<-totalReads-unmappedReads
  counts[["AB_UMIs"]]<-mappedUMIs
  countsTotal[[subsample]]<-unlist(counts)

}

df<-do.call("rbind",countsTotal)
df<-as.data.frame(df)
df$saturation<-1-df$AB_UMIs.CITE_unmapped/df$AB_reads.CITE_unmapped

pdf(paste0(figDir,sampleName,".saturation.pdf"),width=4,height=4)
qplot(1:10,df$saturation,ylim=c(0,1),ylab="saturation")
dev.off()

###### calculate number of working antibodies #######

cite<-citeCells@assays$ADT@counts

igGs<-cite[grepl("IgG\\d",row.names(cite)),]
citeCells@meta.data$cell.size<-log1p(apply(igGs,2,sum))
citeCells@meta.data$log.median.isotype.count<-log1p(apply(igGs,2,median))


pdf(paste0(figDir,sampleName,"_",treatment,"_violin_cellSize.pdf"),width=8,height=6)
VlnPlot(object = citeCells, features="log.median.isotype.count", group.by="garnett_seurat_cluster_call_major_PC_C_res.1.2",pt.size = 0.5)

dev.off()


pdf(paste0(figDir,sampleName,"_",treatment,"_cellSize.pdf"),width=12,height=16)
DimPlot(object = citeCells, reduction="UMAPC",group.by="cell.size")
dev.off()


break()



df<-do.call("rbind",lossTotal)
dfPerAb<-apply(df,2,function(x){x/unlist(abCounts)})
dataM<-melt(dfPerAb)
colnames(dataM)<-c("sample","class","counts")
dataM$class<-gsub("\\..*","",dataM$class)
dataM$class<-factor(dataM$class,levels=names(loss))

pdf(paste0(figDir,"lossAll_counts_perAB.pdf"),width=6,height=6)
p<-ggplot(dataM,aes(sample,counts,fill=class))
p<-p+geom_bar(stat="identity")
p<-p+theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
p<-p+scale_fill_manual(values=cols)

p
dev.off()


df<-do.call("rbind",lossTotal)
dfPerAb<-apply(df,2,function(x){x/unlist(abCounts)})
dataM<-melt(dfPerAb)
colnames(dataM)<-c("sample","class","counts")
dataM$class<-gsub("\\..*","",dataM$class)
dataM$class<-factor(dataM$class,levels=names(loss))

pdf(paste0(figDir,"lossAll_counts_perAB.pdf"),width=6,height=6)
p<-ggplot(dataM,aes(sample,counts,fill=class))
p<-p+geom_bar(stat="identity")
p<-p+theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
p<-p+scale_fill_manual(values=cols)

p
dev.off()

cellNums
df<-do.call("rbind",lossTotal)
dfPerAb<-apply(df,2,function(x){a=x/unlist(abCounts);a=a/unlist(cellNums)})
dataM<-melt(dfPerAb)
colnames(dataM)<-c("sample","class","counts")
dataM$class<-gsub("\\..*","",dataM$class)
dataM$class<-factor(dataM$class,levels=names(loss))

pdf(paste0(figDir,"lossAll_counts_perCell_perAB.pdf"),width=6,height=6)
p<-ggplot(dataM,aes(sample,counts,fill=class))
p<-p+geom_bar(stat="identity")
p<-p+theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
p<-p+scale_fill_manual(values=cols)

p
dev.off()

df<-do.call("rbind",lossTotal)
dfPercent<-t(apply(df,1,function(x){x/sum(x)}))
dataM<-melt(dfPercent)
colnames(dataM)<-c("sample","class","counts")
dataM$class<-gsub("\\..*","",dataM$class)
dataM$class<-factor(dataM$class,levels=names(loss))

pdf(paste0(figDir,"lossAll_counts_percentage.pdf"),width=6,height=6)
p<-ggplot(dataM,aes(sample,counts,fill=class))
p<-p+geom_bar(stat="identity")
p<-p+theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
p<-p+scale_fill_manual(values=cols)

p
dev.off()
























