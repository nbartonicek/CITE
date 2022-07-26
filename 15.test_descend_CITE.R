library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(sctransform)
library(ggplot2)
library(descend)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="MIA"
sampleName="MIA"
treatment="v3"
#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))


matrix_dir = paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_CITE.v3.tsv/umi_count/")
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
cite<-as.data.frame(mat)

rownames(cite) <- paste0("CITE_", rownames(cite))
rownames(cite) <- gsub("-\\w{15}$","",rownames(cite))
cite <-cite[,order(colnames(cite))]
cite <-cite[!grepl("unmapped",row.names(cite)), ]

tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE_vst.Rdata")
load(tSNEFile)
#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells$RNA@data) ]
citeCells[["ADT"]] <- CreateAssayObject(counts = cite)

citeCells <- NormalizeData(citeCells, assay = "ADT", normalization.method = "CLR")
citeCells <- ScaleData(object = citeCells, assay = "ADT")

data=t(cite)
dataM<-melt(data)
colnames(dataM)<-c("cell","AB","count")

pdf(paste0(figDir,"CITE_pool1_test.pdf"),width=12,height=8)
p<-ggplot(dataM,aes(count,group=AB))
p<-p+geom_bar()
p<-p+facet_wrap(~AB,nrow=4)
p
dev.off()

dataM$logCount<-log1p(dataM$count)

pdf(paste0(figDir,"CITE_test_log.pdf"),width=12,height=8)
p<-ggplot(dataM,aes(logCount,group=AB))
p<-p+geom_bar()
p<-p+facet_wrap(~AB,nrow=4)
p
dev.off()


pdf(paste0(figDir,"CITE_",projectname,"_logY_test.pdf"),width=12,height=40)
p<-ggplot(dataM,aes(count,group=AB))
p<-p+geom_bar()
p<-p+facet_wrap(~AB,ncol=6)
p<-p+scale_y_log10()
p<-p+scale_x_log10()
p
dev.off()



for(ab in row.names(cite)){

}



result<-runDescend(cite)

intData<-citeCells@assays$ADT@data
intData<-as.integer(intData)
result<-runDescend(citeCells@assays$ADT@data,scaling.consts=citeCells$nFeature_RNA)

temp<-result

result<-



