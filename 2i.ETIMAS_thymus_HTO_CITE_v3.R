library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
#library(DropletUtils)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="ETIMAS_T21"
sampleName="ETIMAS_T21"
treatment="v3"

dims=40
resolution=0.8
perplexity=20
conditions<-paste(dims,resolution,perplexity,sep="_")



#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))



matrix_dir = paste0(homeDir,"raw_files/",projectname,"/HTO_fastqs/",projectname,"_HTO.v3.tsv/umi_count/")
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



earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_",treatment,"_earlySeurat.Rdata")
if(!file.exists(earlySeuratFile)){
  #dataH <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_feature_bc_matrices/GRCh38")
  dataM <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_feature_bc_matrix")

  #seuratX <- Read10X(dataH)
  seuratM <- Read10X(dataM)
  seuratX<-rbind(seuratM)

  #seuratMat<-as.matrix(seuratX[,colSums(seuratX)>500])
  #this produces 4153 cells (cutoff 500)
  whitelist=read.table(file=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_whitelist.txt"))
  seuratMat=seuratX[,colnames(seuratX) %in% whitelist$V1]
  #this produces 4757 cells  (FDR 0.01)

  #species specific
  #humanIdx<-grepl("GRCh38_",row.names(seuratMat))
  #speciesRatio<-apply(seuratMat,2,function(x){sum(x[humanIdx])/sum(x)})

  #higher proportion of human than mouse
  #pdf(paste0(figDir,"species_ratio_",sampleName,".pdf"),width=12,height=8)
  #plot(density(speciesRatio))
  #dev.off()
  #pdf(paste0(figDir,"species_ratio_hist_",sampleName,".pdf"),width=12,height=8)
  #plot(hist(speciesRatio))
  #dev.off()
  ##eliminate doublets (0 with cutoff of 0.9)
  #seuratMatClean<-seuratMat[,speciesRatio>=0.9|speciesRatio<=0.1]
  seuratMatClean<-seuratMat
  #
  ##species list
  #speciesList<-speciesRatio[speciesRatio>=0.9|speciesRatio<=0.1]
  #speciesNames<-ifelse(speciesList>0.5,"human","mouse")
  #speciesColor<-ifelse(speciesList>0.5,"darkblue","red")
  #save(speciesColor,file="speciesColor.Rdata")
  #save(speciesRatio,file="speciesRatio.Rdata")

  #add the cite data
  seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(hash) ]
  seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(cite) ]


  #
  citeCells <- CreateSeuratObject(
    counts = seuratMatClean, 
    min.cells = 3, 
    min.features = 100, 
    project = sampleName
  )

  save(citeCells,file=earlySeuratFile)
} else {load(earlySeuratFile)}



tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE_vst.Rdata")
if(!file.exists(tSNEFile)){
  pcsData<-as.data.frame(citeCells@reductions$pca@jackstraw@overall.p.values)
  pcs<-pcsData[pcsData$Score<0.05,"PC"]
  citeCells <- FindNeighbors(object = citeCells, dims = pcs, verbose = FALSE)
  citeCells <- FindClusters(citeCells, dims.use = pcs, resolution = 0.8)
  citeCells <- RunTSNE(citeCells, dims.use = pcs)
  save(citeCells,file=tSNEFile)
} else {load(tSNEFile)}



pdf(paste0(figDir,sampleName,"_",treatment,"_tsne.pdf"),width=12,height=8)
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



#genes<-c("nFeature_RNA","percent.mt",paste0(c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
genes<-c("nFeature_RNA","percent.mt",paste0(c("Adam10","Adam17","Bcl2","Bcl2l1","Bcl2l11","Cd28","Cd38","Cd4","Cd44","Cd8","Egr2","Egr3","Entpd1","Fas","Fasl","Gfi1","Gzmb","Hcst","Icos","Ifng","Il15","Il15ra","Il2","Il2rb","Il7r","Itga2","Itgae","Itgal","Klrc1","Klrc2","Klrc3","Klrd1","Klrg1","Klrk1","Nt5e","Prdm1","Prf1","Sell","Tnf")))


pdf(paste0(figDir,sampleName,"_",treatment,"_feature_vst.pdf"),width=12,height=32)
FeaturePlot(object = citeCells, features = genes, ncol = 3, pt.size = 0.5,cols = c("grey", "blue"), 
            reduction = "tsne")

dev.off()
#}


#add CITE data
hash <-hash[,colnames(hash) %in% colnames(citeCells$RNA@data) ]
hash <- hash[1:3,]


positive_quantile=0.999

citeCells[["HTO"]] <- CreateAssayObject(counts = hash)
citeCells <- NormalizeData(object = citeCells, assay = "HTO", normalization.method = "CLR")
citeCells <- HTODemux(citeCells, assay = "HTO", positive.quantile = positive_quantile)
table(citeCells$HTO_classification.global)

Idents(object = citeCells) <- "HTO_maxID"
pdf(paste0(figDir,"ridgeplot_hashv3_",sampleName,".pdf"),width=8,height=6)
RidgePlot(object = citeCells, assay = "HTO", features = rownames(x = citeCells[["HTO"]])[1:3], 
    ncol = 2)
dev.off()

Idents(object = citeCells) <- "HTO_classification"

citeCells <- SetAllIdent(citeCells,"HTO_classification")
cols<-c("darkred", "pink", "pink", 
  "green","lightgreen","darkblue",
                        "gray")

pdf(paste0(figDir,"tsneplot_PQ.",positive_quantile,".",sampleName,".pdf"),width=12,height=8)
TSNEPlot(citeCells,group.by = "HTO_classification", 
         pt.size = 1, 
         label.size = 6,
         cols = cols)
dev.off()

Idents(object = citeCells) <- "HTO_classification.global"
cols<-c("#b2182b", "gray", "#053061")
pdf(paste0(figDir,"hto_classification_global_PQ.",positive_quantile,".",sampleName,".pdf"),width=12,height=8)
TSNEPlot(citeCells,group.by = "HTO_classification.global", 
         pt.size = 1, 
         label.size = 6,
         colors.use = cols)
dev.off()






cite <-cite[,colnames(cite) %in% colnames(citeCells$RNA@data) ]

citeCells[["ADT"]] <- CreateAssayObject(counts = cite)
citeCells <- NormalizeData(citeCells, assay = "ADT", normalization.method = "CLR")
citeCells <- ScaleData(object = citeCells, assay = "ADT")

pdf(paste0(figDir,sampleName,"_",treatment,"_feature_singlet_CITE.pdf"),width=12,height=16)
FeaturePlot(object = citeCells, features = row.names(citeCells$ADT@data), min.cutoff = "q05", max.cutoff = "q90", ncol = 4)
dev.off()

finalObject<-paste0("../project_results/Robjects/",sampleName,"_",treatment,"_CITE_HTO.Rdata")
save(citeCells,file=finalObject)