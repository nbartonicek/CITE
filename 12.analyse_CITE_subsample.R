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


figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))
rObjectsDir=paste0(homeDir,"project_results/Robjects/")
system(paste0("mkdir -p ",rObjectsDir))
annotationDir=paste0(homeDir,"annotation/")
inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530/CITE_fastqs/subsamples"

files=list.files(inDir,pattern="CID4530_CITE",full.names=T)
files=files[grepl("tsv",files)]
#files=files[!grepl("00",files)]


tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE.Rdata")
load(tSNEFile)


countsL<-list()
for(file in files){
  sampleName=gsub(".*CITE_","",basename(file))
  sampleName=gsub(".emptydrop.tsv","",sampleName)
  cat(".")
  #cite <- read.csv(file, sep = ",", header = TRUE, row.names = 1)
  cite <- fread(file)
  cite <- as.data.frame(cite)
  row.names(cite) <- cite[,1]
  cite=cite[,-1]
  rownames(cite) <- paste0("CITE_", rownames(cite))
  cite <- cite[!grepl("bad_struct|no_match|total_reads",row.names(cite)),]
  
  cite <-cite[,order(colnames(cite))]
  if(dim(cite)[1]>0){
    tempDf<-data.frame(cite=names(rowSums(cite)),counts=rowSums(cite))
    countsL[[as.character(sampleName)]]<-tempDf
  }
}

df<-do.call("rbind",countsL)
df$conc=as.numeric(gsub("\\.1\\.C.*|\\.2\\.C.*|\\.3\\.C.*","",row.names(df)))
df$conc[!grepl("CITE",row.names(df))]=as.numeric(gsub("\\.1.*|\\.2.*|\\.3.*","",row.names(df)[!grepl("CITE",row.names(df))]))


nameDF<-data.frame(cite=unique(df$cite),row=rep(1:8,each=4),column=rep(1:4))
nameDF=nameDF[order(nameDF$cite),]
merged<-merge(df,nameDF)
df=merged
pdf("CITE_normal_all.pdf",width=16,height=16)
p<-ggplot(df,aes(conc,counts,colour=cite))
p<-p+geom_smooth()
p<-p+geom_point()
p<-p+scale_x_log10()
p
dev.off()


pdf("CITE_normal_counts.pdf",width=16,height=16)
p<-ggplot(df,aes(conc,counts,colour=cite))
p<-p+geom_smooth(show.legend = FALSE)
p<-p+geom_point(show.legend = FALSE)
p<-p+scale_x_log10()
p<-p+facet_wrap(~cite, ncol=6)+ guides(color=FALSE)
p
dev.off()


############# hashing ############

files=list.files(inDir,pattern="CID4530_HTO",full.names=T)
files=files[grepl("tsv",files)]
files=files[grepl(".1.emptydrop",files)]
files=files[!grepl("00",files)]


#tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE.Rdata")
#load(tSNEFile)


countsL<-list()
countsType<-list()
for(file in files){
  sampleName=gsub(".*HTO_","",basename(file))
  sampleName=gsub(".emptydrop.tsv","",sampleName)
  cat(".")
  #cite <- read.csv(file, sep = ",", header = TRUE, row.names = 1)
  cite <- fread(file)
  cite <- as.data.frame(cite)
  row.names(cite) <- cite[,1]
  cite=cite[,-1]
  #rownames(cite) <- paste0("HTO_", rownames(cite))
  cite <- cite[!grepl("bad_struct|no_match|total_reads",row.names(cite)),]
  
  cite <-cite[,order(colnames(cite))]
  if(dim(cite)[1]>0){
    tempDf<-data.frame(cite=names(rowSums(cite)),counts=rowSums(cite))
    countsL[[as.character(sampleName)]]<-tempDf
  }
  hash=cite
  joint_bcs <- intersect(colnames(citeCells@data),colnames(hash))
  hash <- as.matrix(hash[,joint_bcs])
  hash <- hash[1:2,]

  citeCells <- SetAssayData(citeCells, assay.type = "HTO", slot = "raw.data", new.data = hash)
  citeCells <- NormalizeData(citeCells, assay.type = "HTO", normalization.method = "genesCLR")
  citeCells <- ScaleData(citeCells, assay.type = "HTO", display.progress = FALSE)
  citeCells <- HTODemux(citeCells,assay.type = "HTO",positive_quantile =  0.99,print.output = FALSE)
  print (table(citeCells@meta.data$hto_classification_global))
  countsType[[as.character(sampleName)]]<-table(citeCells@meta.data$hto_classification_global)

}

df<-do.call("rbind",countsL)
df$conc=as.numeric(gsub("\\.1$","",row.names(df)))
df$conc[!grepl("HTO",row.names(df))]=as.numeric(gsub("\\.1.*|\\.2.*|\\.3.*","",row.names(df)[!grepl("HTO",row.names(df))]))
df=df[!grepl("HTO3",row.names(df)),]

#nameDF<-data.frame(cite=unique(df$cite),row=rep(1:8,each=4),column=rep(1:4))
#nameDF=nameDF[order(nameDF$cite),]
#merged<-merge(df,nameDF)
pdf(paste0(figDir,"HASH_normal_all.pdf"),width=8,height=8)
p<-ggplot(df,aes(conc,counts,colour=cite))
p<-p+geom_smooth()
p<-p+geom_point()
#p<-p+scale_x_log10()
p
dev.off()


df<-do.call("rbind",countsType)
df=as.data.frame(df)
df$conc=as.numeric(gsub("\\.1$","",row.names(df)))
dataM<-melt(df)
dataM$conc=unique(df$conc)
dataM<-dataM[1:54,]
pdf(paste0(figDir,"4530_HASH_saturation.pdf"),width=8,height=6)
p<-ggplot(dataM,aes(conc,value,colour=variable))
p<-p+geom_smooth()
p<-p+geom_point()
p<-p+scale_x_log10()+scale_y_log10()
#p<-p+facet_wrap(~variable, ncol=3)
p
dev.off()

########### CITE

files=list.files(inDir,pattern="CID4530_CITE",full.names=T)
files=files[grepl("tsv",files)]
files=files[grepl(".1.emptydrop",files)]
files=files[!grepl("0.00",files)]

igGs<-read.csv(paste0(annotationDir,"igG_content.csv"),header=F)

#tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE.Rdata")
#load(tSNEFile)


SetIfNull <- function(x, default) {
  if(is.null(x = x)){
    return(default)
  } else {
    return(x)
  }
}
MaxN <- function(x, N = 2){
  len <- length(x)
  if (N > len) {
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x, partial = len - N + 1)[len - N + 1]
}
HTODemuxPois<-function (object, assay.type = "HTO", positive_quantile = 0.99, 
    init_centers = NULL, cluster_nstarts = 100, k_function = "clara", 
    nsamples = 100, print.output = TRUE) 
{
    hash_data <- GetAssayData(object = object, assay.type = assay.type)
    hash_raw_data <- GetAssayData(object = object, assay.type = assay.type, 
        slot = "raw.data")[, object@cell.names]
    hash_raw_data <- as.matrix(hash_raw_data)
    ncenters <- SetIfNull(x = init_centers, default = nrow(hash_data) + 
        1)
    if (k_function == "kmeans") {
        hto_init_clusters <- kmeans(x = t(x = GetAssayData(object = object, 
            assay.type = assay.type)), centers = ncenters, nstart = cluster_nstarts)
        object <- SetIdent(object = object, cells.use = names(hto_init_clusters$cluster), 
            ident.use = hto_init_clusters$cluster)
    } else {
        hto_init_clusters <- clara(x = t(x = GetAssayData(object = object, 
            assay.type = assay.type)), k = ncenters, samples = nsamples)
        object <- SetIdent(object = object, cells.use = names(x = hto_init_clusters$clustering), 
            ident.use = hto_init_clusters$clustering)
    }
    object2 <- object
    object2@data <- object2@data[1:10, ]
    hto_averages <- AverageExpression(object = object2, return.seurat = TRUE, 
        show.progress = FALSE)
    average_hto <- GetAssayData(object = hto_averages, assay.type = assay.type, 
        slot = "raw.data")
    hto_discrete <- GetAssayData(object = object, assay.type = assay.type)
    hto_discrete[hto_discrete > 0] <- 0
    for (hto_iter in rownames(x = hash_data)) {
        hto_values <- hash_raw_data[hto_iter, object@cell.names]
        hto_values_use <- hto_values[WhichCells(object = object, 
            ident = which.min(x = average_hto[hto_iter, ]))]
        hto_fit <- fitdist(hto_values_use, method="mme","pois")
        hto_cutoff <- as.numeric(x = quantile(x = hto_fit, probs = positive_quantile)$quantiles[1])
        hto_discrete[hto_iter, names(x = which(x = hto_values > 
            hto_cutoff))] <- 1
        if (print.output) {
            print(paste0("Cutoff for ", hto_iter, " : ", hto_cutoff, 
                " reads"))
        }
    }
    num_hto_positive <- colSums(x = hto_discrete)
    hto_classification_global <- num_hto_positive
    hto_classification_global[num_hto_positive == 0] <- "Negative"
    hto_classification_global[num_hto_positive == 1] <- "Singlet"
    hto_classification_global[num_hto_positive > 1] <- "Doublet"
    donor_id = rownames(x = hash_data)
    hash_max <- apply(X = hash_data, MARGIN = 2, FUN = max)
    hash_maxID <- apply(X = hash_data, MARGIN = 2, FUN = which.max)
    hash_second <- apply(X = hash_data, MARGIN = 2, FUN = MaxN, 
        N = 2)
    hash_maxID <- as.character(x = donor_id[sapply(X = 1:ncol(x = hash_data), 
        FUN = function(x) {
            return(which(x = hash_data[, x] == hash_max[x])[1])
        })])
    hash_secondID <- as.character(x = donor_id[sapply(X = 1:ncol(x = hash_data), 
        FUN = function(x) {
            return(which(x = hash_data[, x] == hash_second[x])[1])
        })])
    hash_margin <- hash_max - hash_second
    doublet_id <- sapply(X = 1:length(x = hash_maxID), function(x) paste(sort(c(hash_maxID[x], 
        hash_secondID[x])), collapse = "_"))
    hto_classification <- hto_classification_global
    hto_classification[hto_classification_global == "Negative"] <- "Negative"
    hto_classification[hto_classification_global == "Singlet"] <- hash_maxID[which(hto_classification_global == 
        "Singlet")]
    hto_classification[hto_classification_global == "Doublet"] <- doublet_id[which(hto_classification_global == 
        "Doublet")]
    classification_metadata <- data.frame(hash_maxID, hash_secondID, 
        hash_margin, hto_classification, hto_classification_global)
    object <- AddMetaData(object = object, metadata = classification_metadata)
    if (print.output) {
        print(x = table(object@meta.data$hto_classification_global))
    }
    object = SetAllIdent(object = object, id = "hto_classification")
    object = SetIdent(object, cells.use = WhichCells(object, 
        subset.name = "hto_classification_global", accept.value = "Doublet"), 
        ident.use = "Doublet")
    object@meta.data$hash_ID = object@ident[rownames(object@meta.data)]
    return(object)
}


sampleName="4530"
treatments=c("normal")
treatment="normal_emptydrop"
if(treatment=="normal_emptydrop"){tag=""}else{tag="-N"}
hash <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/CID",sampleName,tag,"_HTO.tsv"), 
                     sep = ",", header = TRUE, row.names = 1)
rownames(hash) <- gsub("HTO","HTO_", rownames(hash))
hash <- hash[1:(dim(hash)[1]-3),]
hash <-hash[colnames(hash) %in% colnames(citeCells) ]
hash <-hash[,order(colnames(hash))]
    





#getIgG
i=9

file=files[i]
sampleName=gsub(".*CITE_","",basename(file))
sampleName=gsub(".emptydrop.tsv","",sampleName)
cat(".")
#cite <- read.csv(file, sep = ",", header = TRUE, row.names = 1)
cite <- fread(file)
cite <- as.data.frame(cite)
row.names(cite) <- cite[,1]

cite=cite[,-1]
#rownames(cite) <- paste0("HTO_", rownames(cite))
cite <- cite[!grepl("bad_struct|no_match|total_reads",row.names(cite)),]
#rownames(cite) <- paste0("HTO_", 1:dim(cite)[1])
rownames(cite) <- paste0("HTO_", rownames(cite))

cite <-cite[,order(colnames(cite))]
if(dim(cite)[1]>0){
  tempDf<-data.frame(cite=names(rowSums(cite)),counts=rowSums(cite))
  countsL[[as.character(sampleName)]]<-tempDf
}
hash=cite
joint_bcs <- intersect(colnames(citeCells@data),colnames(hash))
hash <- as.matrix(hash[,joint_bcs])
igGHash <-hash[grepl("IgG",row.names(hash)),]

countsL<-list()
countsType<-list()
for (i in 1:length(files)){
  file=files[i]
  sampleName=gsub(".*CITE_","",basename(file))
  sampleName=gsub(".emptydrop.tsv","",sampleName)
  cat(".")
  #cite <- read.csv(file, sep = ",", header = TRUE, row.names = 1)
  cite <- fread(file)
  cite <- as.data.frame(cite)
  row.names(cite) <- cite[,1]

  cite=cite[,-1]
  #rownames(cite) <- paste0("HTO_", rownames(cite))
  cite <- cite[!grepl("bad_struct|no_match|total_reads",row.names(cite)),]
  #rownames(cite) <- paste0("HTO_", 1:dim(cite)[1])
  rownames(cite) <- paste0("HTO_", rownames(cite))

  cite <-cite[,order(colnames(cite))]
  if(dim(cite)[1]>0){
    tempDf<-data.frame(cite=names(rowSums(cite)),counts=rowSums(cite))
    countsL[[as.character(sampleName)]]<-tempDf
  }
  hash=cite
  joint_bcs <- intersect(colnames(citeCells@data),colnames(hash))
  hash <- as.matrix(hash[,joint_bcs])

  res<-list()
  #for(j in 1:dim(hash)[1]){
  indexes<-list()
  for(j in c(19)){
      ab=gsub("HTO_","",row.names(hash))[j]
    if(!grepl("IgG",ab)){
      cat(".")
      igg=as.character(igGs[,2][igGs[,1]==ab])
      hashComp=paste0("HTO_",c(ab,igg))
      tempHash <- as.data.frame(hash[row.names(hash) %in% hashComp,])
      tempHash[hashComp[2],]<-igGHash[hashComp[2],]
      citeCells <- SetAssayData(citeCells, assay.type = "HTO", slot = "raw.data", new.data = tempHash)
      citeCells <- NormalizeData(citeCells, assay.type = "HTO", normalization.method = "genesCLR")
      #citeCells <- NormalizeData(citeCells, assay.type = "HTO")
      citeCells <- ScaleData(citeCells, assay.type = "HTO",display.progress = FALSE)
      citeCells <- HTODemuxPois(citeCells,assay.type = "HTO",k_function = "kmeans",positive_quantile =  0.99,print.output = FALSE)
      print (table(citeCells@meta.data$hto_classification_global))
      res[[j]]<-table(citeCells@meta.data$hto_classification_global)
    }
    indexes[[as.character(j)]]=j
  }
  resDF<-do.call("rbind",res)
  row.names(resDF)<-row.names(hash)[unlist(indexes)]
  countsType[[as.character(sampleName)]]<-resDF
}





df<-do.call("rbind",countsType)
df=as.data.frame(df)
#row.names(df)=paste0(row.names(df),"_",names(countsType))
#df$conc=as.numeric(gsub("\\.1$","",names(countsType)))
dataM<-melt(df)
dataM$conc=rep(as.numeric(gsub("\\.1$","",names(countsType))),each=length(row.names(resDF)))
dataM$cite=gsub("HTO_","",row.names(resDF))

#dataM<-dataM[1:108,]
pdf(paste0(figDir,"4530_CITE_saturation_test.pdf"),width=8,height=6)
p<-ggplot(dataM,aes(conc,value,colour=variable))
p<-p+geom_smooth()
p<-p+geom_point()
p<-p+scale_x_log10()
#+scale_y_log10()
p<-p+facet_wrap(~cite, ncol=1)
p
dev.off()






