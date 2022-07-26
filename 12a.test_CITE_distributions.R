library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(data.table)
library(ggplot2)
library(cluster)
library(fitdistrplus)


homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4530"
sampleName="4530"
treatments=c("normal")
treatment="normal_emptydrop"


figDir=paste0(homeDir,"project_results/figures/CITE_optimisation/")
system(paste0("mkdir -p ",figDir))
rObjectsDir=paste0(homeDir,"project_results/Robjects/")
system(paste0("mkdir -p ",rObjectsDir))
annotationDir=paste0(homeDir,"annotation/")
inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530/CITE_fastqs/subsamples"

files=list.files(inDir,pattern="CID4530_CITE",full.names=T)
files=files[grepl("tsv",files)]
#files=files[!grepl("00",files)]


tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE.Rdata")
#load(tSNEFile)
clusters<-citeCells@ident

########### CITE



files=list.files(inDir,pattern="CID4530_CITE",full.names=T)
files=files[grepl("tsv",files)]
subFiles=files[grepl("0\\.(1|5)\\.1.emptydrop",files)]
igGs<-read.csv(paste0(annotationDir,"igG_content.csv"),header=F)

prepareCite<-function(file){
  cite <- fread(subFile)
  cite <- as.data.frame(cite)
  row.names(cite) <- cite[,1]

  cite=cite[,-1]
  #rownames(cite) <- paste0("HTO_", rownames(cite))
  cite <- cite[!grepl("bad_struct|no_match|total_reads",row.names(cite)),]
  #rownames(cite) <- paste0("HTO_", 1:dim(cite)[1])
  #rownames(citeSubset) <- paste0("HTO_", rownames(citeSubset))
  joint_bcs <- intersect(colnames(citeCells@data),colnames(cite))
  cite <- as.matrix(cite[,joint_bcs])
}

normalizeCite<-function(df){
  res<-list()
  citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = df)
  citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
  citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)
  
  #raw
  res[["raw"]] <- citeCells@assay$CITE@raw.data
  # Normalized
  res[["normalized"]] <- citeCells@assay$CITE@data
  # Scaled
  res[["scaled"]] <- citeCells@assay$CITE@scale.data
  return(res)
}

cites<-list()
oriFile="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530/CITE_fastqs/CID4530_CITE.tsv"
citeCountsL<-list()
for(subFile in c(subFiles,oriFile)){
  if(!grepl("CID4530_CITE.tsv",subFile)){
    sampleName=gsub("CID4530_CITE_","",basename(subFile))
    sampleName=gsub(".1.emptydrop.tsv","",sampleName)
    cat(sampleName)
    cat("\n")
  } else {
    sampleName="1.0"
    cat("1.0\n")
  }
  cite<-prepareCite(subFile)
  citeCountsL[[sampleName]]<-normalizeCite(cite)
}

#for a known antibody plot all diagnostic plots for 1/2/5/10?
#1. raw/normalized/scaled vs IgG
#2. expression

#extract data
abs<-c("CD4","CD45","CD183")

for(ab in abs){

igg=as.character(igGs[,2][igGs[,1]==ab])
citeIDs=paste0("",c(ab,igg))

covData<-list()
for(coverage in names(citeCountsL)){
  abDataL<-list()
  for(type in names(citeCountsL[[1]])){
    dataset=citeCountsL[[coverage]][[type]]
    abDataL[[type]] <- as.data.frame(dataset[row.names(dataset) %in% citeIDs,])
  }
  df<-do.call("rbind",abDataL)
  covData[[coverage]]<-df
}

df<-do.call("rbind",covData)
df$coverage<-rep(names(citeCountsL),each=6)
df$type<-rep(names(citeCountsL[[1]]),each=2)
df$ab<-citeIDs
dataM<-melt(df)

subsample<-unique(dataM$variable)
dataMS<-dataM[dataM$variable %in% subsample,]
dataMS$type=factor(dataMS$type,levels=names(citeCountsL[[1]]))
dataMS$ab=factor(dataMS$ab,levels=rev(citeIDs))


variableOrder<-as.character(dataMS$variable[dataMS$coverage=="1.0"&dataMS$ab==ab&dataMS$type=="raw"])[order(dataMS$value[dataMS$coverage=="1.0"&dataMS$ab==ab&dataMS$type=="raw"],decreasing=T)]
dataMS$variable=factor(dataMS$variable,levels=variableOrder)

#first plot only raw for one AB
dataShort=dataMS[dataMS$coverage=="1.0"&dataMS$type=="raw",]

cellNum=dim(citeCells@data)[2] 

pdf(paste0(figDir,ab,"_raw.pdf"),width=12,height=8)
p<-ggplot(dataShort,aes(variable,value,group=ab,color=ab))
#p<-p+geom_point()
p<-p+geom_line(size=1.2)
#p<-p+scale_x_discrete("variable",breaks=interaction(dataMS$variable,dataMS$ab),labels=dataMS$variable) 
#p<-p+facet_grid(type~coverage, scales = "free")
p<-p+scale_color_manual(values=c("steelblue","darkred"))
p<-p+ ggtitle(paste0("Raw CITE data for ",ab)) 
p<-p+ xlab("cells") 
p<-p+ ylab("counts") 
#p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
#p<-p+scale_x_discrete(name="cells", breaks=500*(1:(1+as.integer(cellNum/500))), labels=500*(1:(1+as.integer(cellNum/500))))
p<-p+scale_x_discrete(labels=500*(0:(as.integer(cellNum/500))),breaks = levels(dataShort$variable)[c(T, rep(F, 499))])
print(p)
dev.off()

#then plot the appropriate cluster on top

df<-do.call("rbind",covData)
df$coverage<-rep(names(citeCountsL),each=6)
df$type<-rep(names(citeCountsL[[1]]),each=2)
df$ab<-citeIDs
dataM<-melt(df)

subsample<-unique(dataM$variable)
dataMS<-dataM[dataM$variable %in% subsample,]
dataMS$type=factor(dataMS$type,levels=names(citeCountsL[[1]]))
dataMS$ab=factor(dataMS$ab,levels=rev(citeIDs))


variableOrder<-as.character(dataMS$variable[dataMS$coverage=="1.0"&dataMS$ab==ab&dataMS$type=="raw"])[order(dataMS$value[dataMS$coverage=="1.0"&dataMS$ab==ab&dataMS$type=="raw"],decreasing=T)]
dataMS$variable=factor(dataMS$variable,levels=variableOrder)

#first plot only raw for one AB
dataShort=dataMS[dataMS$coverage=="1.0"&dataMS$type=="raw",]
isCluster4<-names(citeCells@ident)[citeCells@ident=="4"]

cellNum=dim(citeCells@data)[2] 
dataShort$cluster="other"
dataShort$cluster[(dataShort$variable %in% isCluster4)]="T-cells"
dataShort$cluster[grepl("Ig",dataShort$ab)]="other"

pdf(paste0(figDir,ab,"_raw_annotatedCluster.pdf"),width=12,height=8)
p<-ggplot(dataShort,aes(variable,value,group=ab,color=ab))
p<-p+geom_line(size=1.2)
p<-p+geom_point(aes(colour=cluster),size=0.4)
#p<-p+scale_x_discrete("variable",breaks=interaction(dataMS$variable,dataMS$ab),labels=dataMS$variable) 
#p<-p+facet_grid(type~coverage, scales = "free")
p<-p+scale_color_manual(values=c("darkred","steelblue","NA","orange"))
p<-p+ ggtitle(paste0("Raw CITE data for ",ab)) 
p<-p+ xlab("cells") 
p<-p+ ylab("counts") 
#p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
#p<-p+scale_x_discrete(name="cells", breaks=500*(1:(1+as.integer(cellNum/500))), labels=500*(1:(1+as.integer(cellNum/500))))
p<-p+scale_x_discrete(labels=500*(0:(as.integer(cellNum/500))),breaks = levels(dataShort$variable)[c(T, rep(F, 499))])
print(p)
dev.off()


#plot the averages

pdf(paste0(figDir,ab,"_raw_smooth.pdf"),width=12,height=8)
p<-ggplot(dataShort,aes(variable,value,group=ab,color=ab))
p<-p+geom_smooth(size=1.2)
p<-p+geom_point(aes(colour=cluster),size=1)
#p<-p+scale_x_discrete("variable",breaks=interaction(dataMS$variable,dataMS$ab),labels=dataMS$variable) 
#p<-p+facet_grid(type~coverage, scales = "free")
p<-p+scale_color_manual(values=c("darkred","steelblue","NA","orange"))
p<-p+ ggtitle(paste0("Raw CITE data for ",ab)) 
p<-p+ xlab("cells") 
p<-p+ ylab("counts") 
#p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
#p<-p+scale_x_discrete(name="cells", breaks=500*(1:(1+as.integer(cellNum/500))), labels=500*(1:(1+as.integer(cellNum/500))))
p<-p+scale_x_discrete(labels=500*(0:(as.integer(cellNum/500))),breaks = levels(dataShort$variable)[c(T, rep(F, 499))])
print(p)
dev.off()



dataShort$cluster="other"
dataShort$cluster[!(dataShort$variable %in% isCluster4)]="non-T-cells"
dataShort$cluster[grepl("Ig",dataShort$ab)]="other"

pdf(paste0(figDir,ab,"_raw_annotatedClusterOpposite.pdf"),width=12,height=8)
p<-ggplot(dataShort,aes(variable,value,group=ab,color=ab))
p<-p+geom_line(size=1.2)
p<-p+geom_point(aes(colour=cluster),size=0.4)
#p<-p+scale_x_discrete("variable",breaks=interaction(dataMS$variable,dataMS$ab),labels=dataMS$variable) 
#p<-p+facet_grid(type~coverage, scales = "free")
p<-p+scale_color_manual(values=c("darkred","steelblue","darkblue","NA"))
p<-p+ ggtitle(paste0("Raw CITE data for ",ab)) 
p<-p+ xlab("cells") 
p<-p+ ylab("counts") 
#p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
#p<-p+scale_x_discrete(name="cells", breaks=500*(1:(1+as.integer(cellNum/500))), labels=500*(1:(1+as.integer(cellNum/500))))
p<-p+scale_x_discrete(labels=500*(0:(as.integer(cellNum/500))),breaks = levels(dataShort$variable)[c(T, rep(F, 499))])
print(p)
dev.off()


#raw, normalized, scaled


#first plot only raw for one AB
dataShort=dataMS[dataMS$coverage=="1.0",]
isCluster4<-names(citeCells@ident)[citeCells@ident=="4"]

cellNum=dim(citeCells@data)[2] 
dataShort$cluster="other"
dataShort$cluster[(dataShort$variable %in% isCluster4)]="T-cells"
dataShort$cluster[grepl("Ig",dataShort$ab)]="other"

pdf(paste0(figDir,ab,"_raw_annotatedCluster_datatypes.pdf"),width=12,height=8)
p<-ggplot(dataShort,aes(variable,value,group=ab,color=ab))
p<-p+geom_line(size=1.2)
p<-p+geom_point(aes(colour=cluster),size=0.4)
#p<-p+scale_x_discrete("variable",breaks=interaction(dataMS$variable,dataMS$ab),labels=dataMS$variable) 
#p<-p+facet_grid(type~coverage, scales = "free")
p<-p+scale_color_manual(values=c("darkred","steelblue","NA","orange"))
p<-p+ ggtitle(paste0("CITE data for ",ab)) 
p<-p+ xlab("cells") 
p<-p+ ylab("counts") 
#p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
#p<-p+scale_x_discrete(name="cells", breaks=500*(1:(1+as.integer(cellNum/500))), labels=500*(1:(1+as.integer(cellNum/500))))
p<-p+scale_x_discrete(labels=500*(0:(as.integer(cellNum/500))),breaks = levels(dataShort$variable)[c(T, rep(F, 499))])
p<-p+theme(axis.text.x=element_text(size=8))
p<-p+facet_grid(type~., scales = "free")
print(p)
dev.off()


#plot the averages

pdf(paste0(figDir,ab,"_raw_smooth_datatypes.pdf"),width=12,height=8)
p<-ggplot(dataShort,aes(variable,value,group=ab,color=ab))
p<-p+geom_smooth(size=1.2)
p<-p+geom_point(aes(colour=cluster),size=1)
#p<-p+scale_x_discrete("variable",breaks=interaction(dataMS$variable,dataMS$ab),labels=dataMS$variable) 
#p<-p+facet_grid(type~coverage, scales = "free")
p<-p+scale_color_manual(values=c("darkred","steelblue","NA","orange"))
p<-p+ ggtitle(paste0("CITE data for ",ab)) 
p<-p+ xlab("cells") 
p<-p+ ylab("counts") 
#p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
#p<-p+scale_x_discrete(name="cells", breaks=500*(1:(1+as.integer(cellNum/500))), labels=500*(1:(1+as.integer(cellNum/500))))
p<-p+scale_x_discrete(labels=500*(0:(as.integer(cellNum/500))),breaks = levels(dataShort$variable)[c(T, rep(F, 499))])
p<-p+facet_grid(type~., scales = "free")
p<-p+theme(axis.text.x=element_text(size=8))
print(p)
dev.off()





#first plot only raw for one AB
dataShort=dataMS
isCluster4<-names(citeCells@ident)[citeCells@ident=="4"]

cellNum=dim(citeCells@data)[2] 
dataShort$cluster="other"
dataShort$cluster[(dataShort$variable %in% isCluster4)]="T-cells"
dataShort$cluster[grepl("Ig",dataShort$ab)]="other"

pdf(paste0(figDir,ab,"_raw_annotatedCluster_datatypes_coverage.pdf"),width=12,height=8)
p<-ggplot(dataShort,aes(variable,value,group=ab,color=ab))
p<-p+geom_line(size=1.2)
p<-p+geom_point(aes(colour=cluster),size=0.4)
#p<-p+scale_x_discrete("variable",breaks=interaction(dataMS$variable,dataMS$ab),labels=dataMS$variable) 
#p<-p+facet_grid(type~coverage, scales = "free")
p<-p+scale_color_manual(values=c("darkred","steelblue","NA","orange"))
p<-p+ ggtitle(paste0("CITE data for ",ab)) 
p<-p+ xlab("cells") 
p<-p+ ylab("counts") 
#p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
#p<-p+scale_x_discrete(name="cells", breaks=500*(1:(1+as.integer(cellNum/500))), labels=500*(1:(1+as.integer(cellNum/500))))
p<-p+scale_x_discrete(labels=500*(0:(as.integer(cellNum/500))),breaks = levels(dataShort$variable)[c(T, rep(F, 499))])
p<-p+facet_grid(type~coverage, scales = "free")
p<-p+theme(axis.text.x=element_text(size=8))
print(p)
dev.off()


#plot the averages

pdf(paste0(figDir,ab,"_raw_smooth_datatypes_coverage.pdf"),width=12,height=8)
p<-ggplot(dataShort,aes(variable,value,group=ab,color=ab))
p<-p+geom_smooth(size=1.2)
p<-p+geom_point(aes(colour=cluster),size=1)
#p<-p+scale_x_discrete("variable",breaks=interaction(dataMS$variable,dataMS$ab),labels=dataMS$variable) 
#p<-p+facet_grid(type~coverage, scales = "free")
p<-p+scale_color_manual(values=c("darkred","steelblue","NA","orange"))
p<-p+ ggtitle(paste0("CITE data for ",ab)) 
p<-p+ xlab("cells") 
p<-p+ ylab("counts") 
p<-p+theme(axis.text.x=element_text(size=8))
#p<-p+scale_x_discrete(name="cells", breaks=500*(1:(1+as.integer(cellNum/500))), labels=500*(1:(1+as.integer(cellNum/500))))
p<-p+scale_x_discrete(labels=500*(0:(as.integer(cellNum/500))),breaks = levels(dataShort$variable)[c(T, rep(F, 499))])
p<-p+facet_grid(type~coverage, scales = "free")
print(p)
dev.off()

}
















