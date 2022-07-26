library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
#library(DropletUtils)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4T1_3x"
sampleName="4T1_3x"
treatment="v3"

#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))

citeCells<-readRDS("../raw_files/4T1_3x/outs/03_seurat_object_processed.RData")


matrix_dir = paste0(homeDir,"raw_files/",projectname,"/HTO_fastqs/",projectname,"_HTO.v3/umi_count/")
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
colnames(hash)<-gsub("^","4T1x3_",colnames(hash))


matrix_dir = paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_CITE.v3/umi_count/")
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
colnames(hash)<-gsub("^","4T1x3_",colnames(hash))
cite <-cite[,order(colnames(citeCells$RNA@data))]




pdf(paste0(figDir,sampleName,"_",treatment,"_tsne.pdf"),width=12,height=8)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F","lightgreen","darkblue","lightblue","darkblue","orange","purple","gray","green","black","pink")
DimPlot(object = citeCells, reduction.use = "umap",
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         cols = cols
)
dev.off()



#genes<-c("nFeature_RNA","percent.mt",paste0(c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
genes<-c("nFeature_RNA","percent.mt",paste0(c("Adam10","Adam17","Bcl2","Bcl2l1","Bcl2l11","Cd28","Cd38","Cd4","Cd44","Cd8","Egr2","Egr3","Entpd1","Fas","Fasl","Gfi1","Gzmb","Hcst","Icos","Ifng","Il15","Il15ra","Il2","Il2rb","Il7r","Itga2","Itgae","Itgal","Klrc1","Klrc2","Klrc3","Klrd1","Klrg1","Klrk1","Nt5e","Prdm1","Prf1","Sell","Tnf")))


pdf(paste0(figDir,sampleName,"_",treatment,"_feature_vst.pdf"),width=12,height=24)
FeaturePlot(object = citeCells, reduction = "UMAPA",features = genes, ncol = 3, pt.size = 0.5,cols = c("grey", "blue"))

dev.off()
#}


#add CITE data
hash <-hash[,colnames(hash) %in% colnames(citeCells$RNA@data) ]
hash <-hash[,order(colnames(hash))]
hash <- hash[c(1:3),]


positive_quantile=0.99

citeCells[["HTO"]] <- CreateAssayObject(counts = hash)
citeCells <- NormalizeData(object = citeCells, assay = "HTO", normalization.method = "CLR")
citeCells <- HTODemux(citeCells, assay = "HTO", positive.quantile = positive_quantile)
table(citeCells$HTO_classification.global)

Idents(object = citeCells) <- "HTO_maxID"

nHash<-dim(hash)[1]
nCol=nHash %% 3
nRow=nHash %/% 3
pdf(paste0(figDir,"ridgeplot_hashv3_",sampleName,".pdf"),width=8,height=4*(nRow+1))
RidgePlot(object = citeCells, assay = "HTO", features = rownames(x = citeCells[["HTO"]])[1:nHash], 
    ncol = 1)
dev.off()

Idents(object = citeCells) <- "HTO_classification"

cols<-c("darkred", "pink", "pink", "pink", 
  "green","lightgreen", "lightgreen","darkblue","lightblue","orange",
                        "gray")
cols<-c("#b2182b", "#1F618D", "pink","#053061", 
                        "#1c9099", "#f4a582","gray")

pdf(paste0(figDir,"umap_PQ.",positive_quantile,".",sampleName,".pdf"),width=12,height=8)
DimPlot(citeCells,reduction.use = "umap",group.by = "HTO_classification", 
         pt.size = 1, 
         label.size = 6,
         cols = cols)
dev.off()

Idents(object = citeCells) <- "HTO_classification.global"
cols<-c("#b2182b", "gray", "#053061")
pdf(paste0(figDir,"hto_classification_global_PQ.",positive_quantile,".",sampleName,".pdf"),width=12,height=8)
DimPlot(citeCells,reduction.use = "umap",group.by = "HTO_classification.global", 
         pt.size = 1, 
         label.size = 6,
         colors.use = cols)
dev.off()





colnames(cite)<-gsub("^","4T1x3_",colnames(cite))

cite <-cite[,colnames(cite) %in% colnames(citeCells$RNA@data) ]
cite <-cite[!grepl("unmapped",row.names(cite)), ]
cite <-cite[,order(colnames(cite))]
row.names(cite)<-gsub("_","-",row.names(cite))

citeCells[["ADT"]] <- CreateAssayObject(counts = cite)
citeCells <- NormalizeData(citeCells, assay = "ADT", normalization.method = "CLR")
citeCells <- ScaleData(object = citeCells, assay = "ADT")


pdf(paste0(figDir,"CITE_",sampleName,"_all.pdf"),width=12,height=16)
FeaturePlot(citeCells, features = row.names(cite)[1:22], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, reduction = "UMAPA",cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()



finalObject<-paste0("../project_results/Robjects/",sampleName,"_",treatment,"_CITE_HTO.Rdata")
save(citeCells,file=finalObject)