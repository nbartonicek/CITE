library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(R.utils)
library(ggfortify)

dims=40
resolution=0.8
perplexities=c(5,10,15,20,50)
perplexity=20


args <- R.utils::commandArgs(asValues=TRUE)
if (!is.null(args[["dims"]])){dims = as.numeric(args$dims)} 
if (!is.null(args[["resolution"]])){resolution = as.numeric(args$resolution)} 
if (!is.null(args[["perplexity"]])){perplexity = as.numeric(args$perplexity)} 


conditions<-paste(dims,resolution,perplexity,sep="_")
cat(conditions)

homeDir="/share/ScratchGeneral/nenbar/projects/CITE/"
homeDir="../../../projects/CITE/"

figDir=paste0(homeDir,"project_results/figures/")
projectname="4515"
sampleName="4515_CITE"
treatment="test"
#treatment="mt0.5_nGene7000"
testFigDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")

system(paste0("mkdir -p ",testFigDir))

cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/4515_CITE.tsv"), 
                       sep = ",", header = TRUE, row.names = 1)
rownames(cite) <- paste0("CITE_", rownames(cite))
cite <- cite[1:(dim(cite)[1]-3),]
cite <-cite[,order(colnames(cite))]


tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_tSNE.Rdata")
load(tSNEFile)
#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)

#load information on the IgG 
igGs<-read.csv(paste0(homeDir,"scripts/igG_content.csv"),header=F)

#get norm data
norm_cite=citeCells@assay$CITE@data

#Find correlation coefficient for each antibody and the appropriate IgG
igGL<-split(igGs,igGs$V2)
results<-list()

igGNames=names(igGL)

for(igName in names(igGL)){

  #for each member calculate correlation coefficient vs the right ab, and vs others
  for(citeAb in igGL[[igName]]$V1){
    citeAb=as.character(citeAb)
    corrs=cor(norm_cite[citeAb,],norm_cite[igName,])
    otherIgGs<-igGNames[!igGNames %in% igName]

    for(otherIgG in otherIgGs){
      corrs=c(corrs,cor(norm_cite[citeAb,],norm_cite[otherIgG,]))
    }
    results[[citeAb]]<-corrs
  }
}


#Find correlation coefficient for each antibody and the appropriate IgG
igGL<-split(igGs,igGs$V2)
results<-list()

igGNames=names(igGL)

for(igName in names(igGL)){

  #for each member calculate correlation coefficient vs the right ab, and vs others
  for(citeAb in igGL[[igName]]$V1){
    citeAb=as.character(citeAb)
    corr1=cor(norm_cite[citeAb,],norm_cite[igGNames[1],])
    corr2=cor(norm_cite[citeAb,],norm_cite[igGNames[2],])
    corr3=cor(norm_cite[citeAb,],norm_cite[igGNames[3],])
    corrs<-c(corr1,corr2,corr3)
    names(corrs)=igGNames
    results[[citeAb]]<-corrs
  }
}

df<-do.call("rbind",results)
df=as.data.frame(df)
df$id=row.names(df)

row.names(norm_cite)=gsub("CITE_","",row.names(norm_cite))
pdf(paste0(figDir,"pca_cite.pdf"),width=12,height=8)
autoplot(prcomp(norm_cite),label = TRUE)
dev.off()


pdf(paste0(figDir,"IgG_binding.pdf"),width=8,height=8)
par(mfrow=c(2,2))
plot(norm_cite[igGNames[1],],norm_cite[igGNames[2],])
plot(norm_cite[igGNames[1],],norm_cite[igGNames[3],])
plot(norm_cite["CITE_CTLA-4",],norm_cite[igGNames[2],])
points(norm_cite["CITE_CTLA-4",],norm_cite[igGNames[1],],col="red")
points(norm_cite["CITE_CTLA-4",],norm_cite[igGNames[3],],col="blue")
plot(norm_cite["CITE_CD25",],norm_cite[igGNames[2],])
points(norm_cite["CITE_CD25",],norm_cite[igGNames[1],],col="red")
points(norm_cite["CITE_CD25",],norm_cite[igGNames[3],],col="blue")
dev.off()

