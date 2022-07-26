library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
#library(DropletUtils)
library(UpSetR)
library(ggplot2)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="3921_3963_4066"
sampleName="3921_3963_4066_HTO"
treatment="v3_noncancer"

dims=40
resolution=0.8
perplexity=20
conditions<-paste(dims,resolution,perplexity,sep="_")



#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))



matrix_dir = paste0(homeDir,"raw_files/",projectname,"/HTO_fastqs/",projectname,"_highMem.tsv/umi_count/")
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
hash<-as.data.frame(as.matrix(mat))

rownames(hash) <- gsub("-\\w{15}$","",rownames(hash))
hash <-hash[grep("HTO",row.names(hash)),]



tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE_vst.Rdata")
load(tSNEFile)



pdf(paste0(figDir,sampleName,"_",treatment,"_tsne_vstGAY.pdf"),width=12,height=8)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F","darkblue","orange")
TSNEPlot(object = citeCells, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         cols = cols
)
dev.off()



genes<-c("nFeature_RNA","percent.mt",paste0(c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,sampleName,"_",treatment,"_feature_vst.pdf"),width=12,height=16)
FeaturePlot(object = citeCells, features = genes, ncol = 3, pt.size = 0.5,cols = c("grey", "blue"), 
            reduction = "tsne")

dev.off()
#}


#add CITE data
hash <-hash[,colnames(hash) %in% colnames(citeCells$RNA@data) ]
citeCells[["HTO"]] <- CreateAssayObject(counts = hash)
citeCells <- NormalizeData(object = citeCells, assay = "HTO", normalization.method = "CLR")
citeCells <- HTODemux(citeCells, assay = "HTO", positive.quantile = 0.999)
table(citeCells$HTO_classification.global)

Idents(object = citeCells) <- "HTO_maxID"
pdf(paste0(figDir,"ridgeplot_hashv3_",sampleName,".pdf"),width=12,height=4)
RidgePlot(object = citeCells, assay = "HTO", features = rownames(x = citeCells[["HTO"]])[1:3], 
    ncol = 3)
dev.off()

Idents(object = citeCells) <- "HTO_classification.global"
citeCells.singlet <- subset(x = citeCells, idents = "Singlet")

Idents(object = citeCells.singlet) <- "HTO_maxID"

#load genotype data


projectname="3921_3963_4066_Blood"
inFolders<-list.files(paste0("../genotyping/",projectname),full.names=T)
inFolders<-inFolders[!grepl("step",inFolders)]

results<-list()
for(inFolder in inFolders){
  files<-list.files(inFolder,full.names=T)
  files<-files[grep("best",files)]
  data=read.table(files,header=T)

  df<-data.frame(bc=data$BARCODE,id=data$BEST)
  df$bc<-gsub("-1","",df$bc)
  df$id<-gsub("5504844386755030621095_UKB_B09_YAN6909A20_","",df$id)
  df$id<-gsub("5504844386755030621095_UKB_B10_YAN6909A21_","",df$id)
  df$id<-gsub("5504844386755030621095_UKB_B11_YAN6909A22_","",df$id)
  df$id<-gsub("_Blood.CEL","",df$id)
  df$id<-gsub("AMB-.*","AMB",df$id)
  df$id<-gsub("SNG-","",df$id)
  df$id<-gsub("DBL-.*","DBL",df$id)

  citeCells.subset <- subset(x = citeCells.singlet, cells=df$bc )
  dfSubset<-df[df$bc %in% colnames(citeCells.subset),]



  hashes<-c("3921","3963","4066")
  names(hashes)<-c("HTO1","HTO2","HTO3")
  res<-list()
  for(i in 1:3){
    bcOl<-dfSubset$bc[dfSubset$id %in% as.character(hashes[i])]

    hto<-citeCells.subset@meta.data$HTO_maxID
    names(hto)<-colnames(citeCells.subset)

    htoOl<-names(hto)[hto %in% names(hashes)[i]]
    res[[i]]<-sum(bcOl %in% htoOl)/length(bcOl)
    cat(sum(bcOl %in% htoOl)/length(bcOl))
    cat("\n")
  }
}


data=read.table("../genotyping/3921_3963_4066_Blood/individual/OTV.best",header=T)


df<-data.frame(bc=data$BARCODE,id=data$BEST)
df$bc<-gsub("-1","",df$bc)
df$id<-gsub("5504844386755030621095_UKB_B09_YAN6909A20_","",df$id)
df$id<-gsub("5504844386755030621095_UKB_B10_YAN6909A21_","",df$id)
df$id<-gsub("5504844386755030621095_UKB_B11_YAN6909A22_","",df$id)
df$id<-gsub("_Blood.CEL","",df$id)
df$id<-gsub("AMB-.*","AMB",df$id)
df$id<-gsub("SNG-","",df$id)
df$id<-gsub("DBL-.*","DBL",df$id)

citeCells.subset <- subset(x = citeCells.singlet, cells=df$bc )
dfSubset<-df[df$bc %in% colnames(citeCells.subset),]

hashes<-c("3921","3963","4066")
names(hashes)<-c("HTO1","HTO2","HTO3")
for(i in 1:3){
  bcOl<-dfSubset$bc[dfSubset$id %in% as.character(hashes[i])]

  hto<-citeCells.subset@meta.data$HTO_maxID
  names(hto)<-colnames(citeCells.subset)

  htoOl<-names(hto)[hto %in% names(hashes)[i]]

  cat(sum(bcOl %in% htoOl)/length(bcOl))
  cat("\n")
}

#Ok these were percentages of overlap, direct one

#How does that look on a tSNE 
#1. What percentage of AMB is AMB.

dfOL<-df[df$bc %in% colnames(citeCells),]
dfOL$class<-sapply(dfOL$id,function(x){ifelse(x=="AMB","Doublet","Singlet")})

citeCells.subset <- subset(x = citeCells, cells=dfOL$bc )
Idents(object = citeCells.subset) <- "HTO_classification.global"


hashTypesL=split(names(Idents(object = citeCells.subset)),Idents(object = citeCells.subset))
genotypeTypesL=split(dfOL$bc,dfOL$class)

overlaps<-list()
for(hashType in names(hashTypesL)){
  overlaps[[hashType]]<-sapply(genotypeTypesL,function(x){sum(hashTypesL[[hashType]] %in% x)})
}



#2. What percentage is identical
#3. Are those that are not more pronounced in certain types of cells
#4. Are all samples the same


#5. What parameters of Hashing and Demuxing can we use to improve this
#does it matter if we use different cutoffs to call singlets in hashing
#does it matter if we use different genotyping parameters, such as a) genotyping with other samples b) using only transcriptome SNPs
#does it matter if I clean with plink
#how does non-blood compare?




#load genotype data
#data=read.table("../genotyping/3921_3963_4066_Blood/step2/3921_3963_4066_Blood.single",header=T)
#data=read.table("../genotyping/3921_3963_4066_Blood/step2/3921_3963_4066_Blood.single",header=T)

#number of used SNPs
annotation<-import("/share/ClusterShare/biodata/contrib/nenbar/genomes/star/hg38_gencode.v27.polyA/genes_merged.gtf")
vcf<-read.table("../genotyping/3921_3963_4066_Blood/step2/3921_3963_4066_Blood.hg38.vcf",header=F)
vcf<-read.table("../genotyping/4513_4530_Blood/step2/4513_4530_Blood.hg38.vcf",header=F)
vcfUnique<-apply(vcf[,10:12],1,function(x){vcf[1,1]==vcf[1,2]&vcf[1,1]==vcf[1,3]})
vcf<-vcf[!vcfUnique,]

grVCF<-GRanges(seqnames=vcf$V1,IRanges(vcf$V2,width=1),strand="*")
mat<-findOverlaps(grVCF,annotation)
length(unique(queryHits(mat)))

pdf(paste0(figDir,"prob_singlet.pdf"),width=6,height=4)
qplot(data$PRB.SNG1,geom="histogram",xlab="Singlet probability",ylab="count")
dev.off()








