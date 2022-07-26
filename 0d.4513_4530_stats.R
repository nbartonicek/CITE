
library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(UpSetR)
library(reshape2)
library(ggplot2)
library(DropletUtils)

homeDir="../../../projects/CITE/"

RobjectsDir=paste0(homeDir,"project_results/Robjects/")

figDir=paste0(homeDir,"project_results/figures/genotyping")
projectname="4513_4530_Blood"
system(paste0("mkdir -p ",figDir))


#calculate percentage overlap with the real thing
#load in barcodes from samples
seurat4513 <- Read10X("../raw_files/4513_4530/outs/4513/")
seurat4530 <- Read10X("../raw_files/4513_4530/outs/4530/")


citeCells4513 <- CreateSeuratObject(
    counts = seurat4513, 
    min.cells = 3, 
    min.features = 100, 
    project = projectname
)
citeCells4530 <- CreateSeuratObject(
    counts = seurat4530, 
    min.cells = 3, 
    min.features = 100, 
    project = projectname
)

bc4513<-colnames(citeCells4513)
bc4530<-colnames(citeCells4530)
bcOL<-bc4513[bc4513 %in% bc4530]

bc<-unique(c(bc4513,bc4530))
bc<-bc[!(bc %in% bcOL)]

bc4513<-bc4513[!(bc4513 %in% bcOL)]
bc4530<-bc4530[!(bc4530 %in% bcOL)]


bc<-paste0(bc,"-1")
write(bc,file="../raw_files/4513_4530/outs/barcodes.tsv")
write(bc4513,file="../raw_files/4513_4530/outs/barcodes4513.tsv")
write(bc4530,file="../raw_files/4513_4530/outs/barcodes4530.tsv")

projectname="4513_4530_Blood"
inFolders<-list.files(paste0("../genotyping/",projectname),full.names=T)
inFolders<-inFolders[grepl("doublet-prior",inFolders)]

results<-list()
cellCount<-list()
for(inFolder in inFolders){
  cat(basename(inFolder))
  cat("\n")
  files<-list.files(inFolder,full.names=T)
  files<-files[grep("best",files)]
  if(length(files)>1){files=files[length(files)]}
  data=read.table(files,header=T)

  df<-data.frame(bc=data$BARCODE,id=data$BEST)
  #df$bc<-gsub("-1","",df$bc)
  df$id<-gsub("5504844386755030621095_UKB_E04_YAN6909A51_","",df$id)
  df$id<-gsub("5504844386755030621095_UKB_B08_YAN6909A19_","",df$id)
  df$id<-gsub("_Blood.*","",df$id)
  df$id<-gsub("AMB-.*","AMB",df$id)
  df$id<-gsub("SNG-","",df$id)
  df$id<-gsub("DBL-.*","DBL",df$id)

  df=df[df$bc %in% totalBC,]

  #calculate which percentage of total cells from df are properly assigned in bc4513/bc4530
  
  count4513<-sum(df$bc[df$id %in% "4513"] %in% bc4513U)
  percent4513<-count4513/length(df$bc[df$id %in% "4513"])

  count4530<-sum(df$bc[df$id %in% "4530"] %in% bc4530U)
  percent4530<-count4530/length(df$bc[df$id %in% "4530"])

  results[[basename(inFolder)]]<-c(percent4513,percent4530)
  cellCount[[basename(inFolder)]]<-c(mean(data$RD.UNIQ[-1]),mean(data$N.SNP[-1]))
}

df<-do.call("cbind",results)
row.names(df)<-c("4513","4530")

dataM<-melt(df)
colnames(dataM)<-c("HTO","type","percentage_overlap")
dataM$HTO=as.character(dataM$HTO)
dataM$type=factor(dataM$type,levels=c("recom","OTV","individual","GWCCG","imputed","imputed_overlap"))

cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F","darkblue","orange")

pdf(paste0(figDir,"/",projectname,"_demux_hash_percent_overlap.pdf"),width=8,height=6, useDingbats=FALSE)
p<-ggplot(dataM,aes(type,percentage_overlap),group=HTO)
p<-p+geom_point(aes(color=HTO),size=5)
p<-p+theme_bw()
p<-p+ylim(0,1)
p
dev.off()




df<-do.call("cbind",cellCount)
row.names(df)<-c("unique_reads","SNPs")

dataM<-melt(df)
colnames(dataM)<-c("class","type","count")
#dataM$class=c("unique_reads","SNPs")
dataM$type=factor(dataM$type,levels=c("recom","OTV","individual","GWCCG","imputed","imputed_overlap"))

pdf(paste0(figDir,"/",projectname,"_demux_hash_read_SNP_counts.pdf"),width=8,height=6, useDingbats=FALSE)
p<-ggplot(dataM,aes(type,count),group=class)
p<-p+geom_point(aes(color=class),size=5)
p<-p+theme_bw()
p<-p+scale_color_manual(values=c("darkred","darkblue"))
p
dev.off()



########### split cells by cluster, calculate average percentage correct
inFolder=inFolders[3]
cat(basename(inFolder))
cat("\n")
files<-list.files(inFolder,full.names=T)
files<-files[grep("best",files)]
if(length(files)>1){files=files[length(files)]}
data=read.table(files,header=T)

df<-data.frame(bc=data$BARCODE,id=data$BEST)
df$id<-gsub("5504844386755030621095_UKB_E04_YAN6909A51_","",df$id)
df$id<-gsub("5504844386755030621095_UKB_B08_YAN6909A19_","",df$id)
df$id<-gsub("_Blood.*","",df$id)
df$id<-gsub("AMB-.*","AMB",df$id)
df$id<-gsub("SNG-","",df$id)
df$id<-gsub("DBL-.*","DBL",df$id)

results<-list()
clusters<-unique(citeCells.subset@meta.data$seurat_clusters)
for(cluster in clusters){
  cells<-colnames(citeCells.subset)[citeCells.subset@meta.data$seurat_clusters==cluster]
  citeCells.subset.cluster <- subset(x = citeCells.singlet, cells=cells )
  dfSubset<-df[df$bc %in% colnames(citeCells.subset.cluster),]
  res<-list()
  for(i in 1:3){
    bcOl<-dfSubset$bc[dfSubset$id %in% as.character(hashes[i])]

    hto<-citeCells.subset.cluster@meta.data$HTO_maxID
    names(hto)<-colnames(citeCells.subset.cluster)

    htoOl<-names(hto)[hto %in% names(hashes)[i]]
    res[[i]]<-sum(bcOl %in% htoOl)/length(bcOl)
    cat(sum(bcOl %in% htoOl)/length(bcOl))
    cat("\n")
  }
  results[[as.character(cluster)]]<-unlist(res)
}
df<-do.call("cbind",results)


row.names(df)<-hashes

dataM<-melt(df)
colnames(dataM)<-c("HTO","cluster","percentage_overlap")
dataM$HTO=as.character(dataM$HTO)
dataM$cluster=as.character(dataM$cluster)

dataM[is.na(dataM)]=0
#dataM$class=c("unique_reads","SNPs")

pdf(paste0(figDir,"/",projectname,"_demux_hash_percentage_per_cluster.pdf"),width=8,height=6, useDingbats=FALSE)
p<-ggplot(dataM,aes(cluster,percentage_overlap),group=HTO)
p<-p+geom_point(aes(color=cluster),size=5)
p<-p+theme_bw()
p<-p+scale_color_manual(values=cols)
p
dev.off()






projectname="4513_4530_Blood"
inFolder=inFolders[1]

results<-list()
cellCount<-list()
precision<-list()
files<-list.files(inFolder,pattern="overlap",full.names=T)
files<-files[grep("best",files)]
files<-files[!grepl("RG",files)]

for(fileName in files){
  cat(basename(fileName))
  cat("\n")
  
  data=read.table(fileName,header=T)

  df<-data.frame(bc=data$BARCODE,id=data$BEST)
  df$id<-gsub("5504844386755030621095_UKB_E04_YAN6909A51_","",df$id)
  df$id<-gsub("5504844386755030621095_UKB_B08_YAN6909A19_","",df$id)
  df$id<-gsub("_Blood.*","",df$id)
  df$id<-gsub("AMB-.*","AMB",df$id)
  df$id<-gsub("SNG-","",df$id)
  df$id<-gsub("DBL-.*","DBL",df$id)

  df=df[df$bc %in% totalBC,]
  count4513<-sum(df$bc[df$id %in% "4513"] %in% bc4513U)
  percent4513<-count4513/length(df$bc[df$id %in% "4513"])

  count4530<-sum(df$bc[df$id %in% "4530"] %in% bc4530U)
  percent4530<-count4530/length(df$bc[df$id %in% "4530"])
  cat(count4513)
  cat("\t")
  cat(count4530)
  cat("\n")  
  results[[basename(fileName)]]<-c(percent4513,percent4530)
  cellCount[[basename(fileName)]]<-c(mean(data$RD.UNIQ[-1]),mean(data$N.SNP[-1]))
  precision[[basename(fileName)]]<-c(count4513+count4530)/dim(df)[1]

}


df<-do.call("cbind",results)
row.names(df)<-c("4513","4530")

dataM<-melt(df)
colnames(dataM)<-c("HTO","subset","percentage_overlap")
dataM$HTO=as.character(dataM$HTO)
dataM$subset=gsub("imputed_overlap_","",dataM$subset)
dataM$subset=gsub(".best","",dataM$subset)
dataM$subset=as.numeric(dataM$subset)/10


pdf(paste0(figDir,"/",projectname,"_demux_hash_subsample__percent_overlap.pdf"),width=8,height=6, useDingbats=FALSE)
p<-ggplot(dataM,aes(subset,percentage_overlap),group=HTO)
p<-p+geom_point(aes(color=HTO),size=5)
p<-p+theme_bw()
p<-p+ylim(0,1)
p<-p+scale_color_manual(values=cols)
p<-p+xlab("Percent subsample")
p
dev.off()


#load in the total sample
if(!file.exists(paste0(RobjectsDir,"4513_4530_miniatlas.Rdata"))){
  finalObject<-paste0(RobjectsDir,"03_seurat_CCA_aligned_processed.Rdata")
  citeCells<-readRDS(finalObject)
  #select 4513/4530
  cells<-colnames(citeCells)[grepl("4513|4530",colnames(citeCells))]
  citeCellsShort<-subset(citeCells,cells=cells)
  citeCells<-citeCellsShort
  saveRDS(citeCells,file=paste0(RobjectsDir,"4513_4530_miniatlas.Rdata"))
} else {citeCells=readRDS(paste0(RobjectsDir,"4513_4530_miniatlas.Rdata"))}


#select cell types and calculate percentage for each cell type 

########### split cells by cluster, calculate average percentage correct
projectname="4513_4530_Blood"
inFolder=inFolders[3]

resultsPercentage<-list()
resultsCounts<-list()
readSNPCount<-list()
precision<-list()
accuracy<-list()


resultsPercentageCluster<-list()
resultsCountsCluster<-list()
readSNPCountCluster<-list()
accuracyCluster<-list()
precisionCluster<-list()

files<-list.files(inFolder,pattern="overlap_",full.names=T)
files<-files[grep("best",files)]
files<-files[!grepl("RG",files)]
clusters<-as.character(unique(citeCells@meta.data$int_cluster_subset_PC_C_res.1.2))

citeCells_BCs<-gsub(".*_","",colnames(citeCells))
bc4513<-citeCells_BCs[grepl("4513",colnames(citeCells))]
bc4530<-citeCells_BCs[grepl("4530",colnames(citeCells))]
bc4513U<-bc4513[!bc4513 %in% bc4530]
bc4530U<-bc4530[!bc4530 %in% bc4513]
totalBC<-c(bc4513U,bc4530U)




for(fileName in files){
  cat(basename(fileName))
  cat("\n")
  
  data=read.table(fileName,header=T)

  df<-data.frame(bc=data$BARCODE,id=data$BEST)
  df$id<-gsub("5504844386755030621095_UKB_E04_YAN6909A51_","",df$id)
  df$id<-gsub("5504844386755030621095_UKB_B08_YAN6909A19_","",df$id)
  df$id<-gsub("_Blood.*","",df$id)
  df$id<-gsub("AMB-.*","AMB",df$id)
  df$id<-gsub("SNG-","",df$id)
  df$id<-gsub("DBL-.*","DBL",df$id)
  df$bc<-gsub("-1","",df$bc)



  df=df[df$bc %in% totalBC,]
  count4513<-sum(df$bc[df$id %in% "4513"] %in% bc4513U)
  count4513_N<-sum(df$bc[df$id %in% "4513"] %in% bc4530U)

  percent4513<-count4513/dim(df)[1]
  percent4513_N<-count4513_N/dim(df)[1]

  count4530<-sum(df$bc[df$id %in% "4530"] %in% bc4530U)
  count4530_N<-sum(df$bc[df$id %in% "4530"] %in% bc4513U)

  percent4530<-count4530/dim(df)[1]
  percent4530_N<-count4530_N/dim(df)[1]

  countDbl<-length(df$bc[df$id %in% c("AMB","DBL")])
  percentDbl<-countDbl/dim(df)[1]

  resultsPercentage[[basename(fileName)]]<-c(percent4513+percent4530,percent4513_N+percent4530_N,percentDbl)
  resultsCounts[[basename(fileName)]]<-c(count4513,count4513_N,count4530,count4530_N,countDbl)

  readSNPCount[[basename(fileName)]]<-c(mean(data$RD.UNIQ[-1]),mean(data$N.SNP[-1]))
  accuracy[[basename(fileName)]]<-c(count4513+count4530)/dim(df)[1]
  precision[[basename(fileName)]]<-c(count4513+count4530)/c(sum(count4513,count4530,count4513_N,count4530_N))




  for(cluster in clusters){
    #cells<-colnames(citeCells)[citeCells@meta.data$seurat_clusters==cluster]
    citeCells_BCsFullname<-colnames(citeCells)[citeCells@meta.data$int_cluster_subset_PC_C_res.1.2==cluster]

    citeCells_BCs<-gsub(".*_","",colnames(citeCells))[citeCells@meta.data$int_cluster_subset_PC_C_res.1.2==cluster]
    bc4513<-citeCells_BCs[grepl("4513",citeCells_BCsFullname)]
    bc4530<-citeCells_BCs[grepl("4530",citeCells_BCsFullname)]
    bc4513U_Clust<-bc4513[!bc4513 %in% bc4530]
    bc4530U_Clust<-bc4530[!bc4530 %in% bc4513]
    totalBCCluster<-c(bc4513U_Clust,bc4530U_Clust)

    dfTemp=df[df$bc %in% totalBCCluster,]
    count4513<-sum(dfTemp$bc[dfTemp$id %in% "4513"] %in% bc4513U_Clust)
    count4513_N<-sum(dfTemp$bc[dfTemp$id %in% "4513"] %in% bc4530U_Clust)

    percent4513<-count4513/dim(dfTemp)[1]
    #percent4513_N<-count4513_N/length(dfTemp$bc[dfTemp$id %in% "4513"])
    #percent4513_N<-count4513_N/length(dfTemp$bc[dfTemp$id %in% "4513"])

    count4530<-sum(dfTemp$bc[dfTemp$id %in% "4530"] %in% bc4530U_Clust)
    count4530_N<-sum(dfTemp$bc[dfTemp$id %in% "4530"] %in% bc4513U_Clust)

    percent4530<-count4530/dim(dfTemp)[1]
    #percent4530_N<-count4530_N/length(dfTemp$bc[dfTemp$id %in% "4530"])
    percent4530_N<-count4530_N/dim(dfTemp)[1]

    countDbl<-length(dfTemp$bc[dfTemp$id %in% c("AMB","DBL")])
    percentDbl<-countDbl/dim(dfTemp)[1]

    num<-gsub("imputed_overlap_","",basename(fileName))
    num<-gsub(".best","",num)

    sampleName=paste0(cluster,"_",num)
    resultsPercentageCluster[[sampleName]]<-c(percent4513+percent4530,percent4513_N+percent4530_N,percentDbl)
    resultsCountsCluster[[sampleName]]<-c(count4513,count4513_N,count4530,count4530_N,countDbl)

    readSNPCountCluster[[sampleName]]<-c(mean(data$RD.UNIQ[-1]),mean(data$N.SNP[-1]))
    accuracyCluster[[sampleName]]<-c(count4513+count4530)/dim(dfTemp)[1]
    precisionCluster[[sampleName]]<-c(count4513+count4530)/c(sum(count4513,count4530,count4513_N,count4530_N))
  }
}

df<-do.call("cbind",results)


#plot of precision/accuracy per UMI count

#take umi
#make true/false
#plot box plots 
#plot dot plots


#make plots of counts 
df<-do.call("cbind",resultsCounts)
row.names(df)<-c("4513_TP","4513_FP","4530_TP","4530_FP","DBL")
dataM<-melt(df)
colnames(dataM)<-c("type","subsample","count")

dataM$subsample<-gsub("imputed_overlap_","",dataM$subsample)
dataM$subsample<-as.numeric(gsub(".best","",dataM$subsample))/10



cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F")

pdf(paste0(figDir,"/",projectname,"_demux_counts_subsample_total.pdf"),width=8,height=6, useDingbats=FALSE)
p<-ggplot(dataM,aes(subsample,count),group=type)
p<-p+geom_point(aes(color=type),size=4)
p<-p+theme_bw()
p<-p+scale_color_manual(values=cols)
p
dev.off()



#make plots of percentage

#make plots of counts 
df<-do.call("cbind",resultsPercentage)
row.names(df)<-c("TP","FP","DBL")
dataM<-melt(df)
colnames(dataM)<-c("type","subsample","count")

dataM$subsample<-gsub("imputed_overlap_","",dataM$subsample)
dataM$subsample<-as.numeric(gsub(".best","",dataM$subsample))/10



cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F")

pdf(paste0(figDir,"/",projectname,"_demux_percentages_subsample_total.pdf"),width=8,height=6, useDingbats=FALSE)
p<-ggplot(dataM,aes(subsample,count),group=type)
p<-p+geom_col(aes(fill=type),position = position_stack())
p<-p+theme_bw()
p<-p+scale_color_manual(values=cols)
p
dev.off()


#make plots of precision
#make plots of counts 
df<-do.call("cbind",accuracy)
row.names(df)<-c("precision")
dataM<-melt(df)
colnames(dataM)<-c("type","subsample","precision")

dataM$subsample<-gsub("imputed_overlap_","",dataM$subsample)
dataM$subsample<-as.numeric(gsub(".best","",dataM$subsample))/10



cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F")

pdf(paste0(figDir,"/",projectname,"_demux_precision_subsample_total.pdf"),width=8,height=6, useDingbats=FALSE)
p<-ggplot(dataM,aes(subsample,precision),group=type)
p<-p+geom_line()
p<-p+ylim(0,1)
p<-p+theme_bw()
p
dev.off()



#########################the same but for clusters

#make plots of counts 
df<-do.call("cbind",resultsCountsCluster)
row.names(df)<-c("4513_TP","4513_FP","4530_TP","4530_FP","DBL")
dataM<-melt(df)
colnames(dataM)<-c("type","cluster","count")
dataM$subsample<-gsub(".*_","",dataM$cluster)

dataM$subsample<-gsub("imputed_overlap_","",dataM$subsample)
dataM$subsample<-as.numeric(gsub(".best","",dataM$subsample))/10

dataM$cluster<-gsub("_\\d.*","",dataM$cluster)


cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F")

pdf(paste0(figDir,"/",projectname,"_demux_counts_subsample_cluster.pdf"),width=10,height=8, useDingbats=FALSE)
p<-ggplot(dataM,aes(subsample,count),group=type)
p<-p+geom_point(aes(color=type),size=2)
p<-p+facet_wrap(~cluster,nrow=4, scales = "free")
p<-p+theme_bw()
p<-p+scale_color_manual(values=cols)
p
dev.off()




#make plots of percentage

#make plots of counts 
df<-do.call("cbind",resultsPercentageCluster)
row.names(df)<-c("TP","FP","DBL")
dataM<-melt(df)
colnames(dataM)<-c("type","cluster","count")
dataM$subsample<-gsub(".*_","",dataM$cluster)


dataM$subsample<-gsub("imputed_overlap_","",dataM$subsample)
dataM$subsample<-as.numeric(gsub(".best","",dataM$subsample))/10

dataM$cluster<-gsub("_\\d.*","",dataM$cluster)


cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F")

pdf(paste0(figDir,"/",projectname,"_demux_percentages_subsample_cluster.pdf"),width=8,height=6, useDingbats=FALSE)
p<-ggplot(dataM,aes(subsample,count),group=type)
p<-p+geom_col(aes(fill=type),position = position_stack())
p<-p+theme_bw()
p<-p+facet_wrap(~cluster,nrow=4, scales = "free")
p<-p+scale_color_manual(values=cols)
p
dev.off()







