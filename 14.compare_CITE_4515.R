library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4515"
sampleName="4515_CITE"
treatment="original"
#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))

cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/4515_CITE.tsv"), 
                       sep = ",", header = TRUE, row.names = 1)
rownames(cite) <- paste0("CITE_", rownames(cite))
cite <- cite[1:(dim(cite)[1]-3),]
citeOld <-cite[,order(colnames(cite))]


homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4515"
sampleName="4515_CITE"
treatment="v3"
#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))


matrix_dir = paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/4515_CITE.v3.tsv/umi_count/")
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "features.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
cite<-as.data.frame(mat)

rownames(cite) <- paste0("CITE_", rownames(cite))
rownames(cite) <- gsub("-\\w{15}$","",rownames(cite))
citeNew <-cite[,order(colnames(cite))]

citeNew<-citeNew[,colnames(citeNew) %in% colnames(citeOld)]
citeOld<-citeOld[,colnames(citeOld) %in% colnames(citeNew)]

pdf(paste0(figDir,"feature_plot_cite_",sampleName,".pdf"),width=12,height=20)
FeaturePlot(citeCells, features.plot = row.names(cite), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()
