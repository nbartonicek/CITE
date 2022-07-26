library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(sctransform)
library(ggplot2)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4515"
sampleName="4515_CITE"
treatment="v3"
#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/forTony/")
system(paste0("mkdir -p ",figDir))


matrix_dir = paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/4515_CITE.v3.tsv/read_count/")
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


earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_earlySeurat3.Rdata")
if(!file.exists(earlySeuratFile)){
  dataH <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/GRCh38")
  whitelist=read.table(file=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_whitelist.txt"))

  seuratH <- Read10X(dataH)
  seuratX <- seuratH
  seuratMat=seuratX[,colnames(seuratX) %in% whitelist$V1]
  #this produces 2646 cells (cutoff 1000)

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
  seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(cite) ]


  #combine the mouse and human
  #combined<-CollapseSpeciesExpressionMatrix(seuratMatClean,prefix.1="GRCh38_",prefix.controls="mm10_",features.controls.toKeep=500)

  citeCells <- CreateSeuratObject(
    counts = seuratMatClean, 
    min.cells = 3, 
    min.features = 100, 
    project = sampleName
  )

  save(citeCells,file=earlySeuratFile)
} else {load(earlySeuratFile)}


citeCells <- PercentageFeatureSet(object = citeCells, pattern = "^MT-", col.name = "percent.mt")


pdf(paste0(figDir,"vlnplotS3_",sampleName,".pdf"),width=12,height=8)
VlnPlot(object = citeCells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

citeCells <- subset(x = citeCells, subset = nFeature_RNA > 250 & percent.mt < 25)
temp<-citeCells
sctranscormFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_sctransform.Rdata")
if(!file.exists(sctranscormFile)){
  citeCells <- SCTransform(object = citeCells, verbose = FALSE)
  citeCells <- RunPCA(object = citeCells, verbose = FALSE,npcs = 100)
  citeCells <- JackStraw(object = citeCells, dims=100,num.replicate=100,verbose=T)
  citeCells <- ScoreJackStraw(citeCells,dims=1:100)
  save(citeCells,file=sctranscormFile)
} else {load(sctranscormFile)}



v3File=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_v3.Rdata")
if(!file.exists(v3File)){
  citeCells <- NormalizeData(object = temp, verbose = FALSE)
  citeCells <- FindVariableFeatures(object = citeCells, selection.method = "dispersion", nfeatures = 2000)
  citeCells <- ScaleData(object = citeCells, verbose = FALSE,vars.to.regress = c("nFeature_RNA","nCount_RNA"))
  citeCells <- RunPCA(object = citeCells, npcs = 100, verbose = FALSE)
  citeCells <- JackStraw(object = citeCells, dims=100,num.replicate=100,verbose=T)
  citeCells <- ScoreJackStraw(citeCells,dims=1:100)

  save(citeCells,file=v3File)
} else {load(v3File)}





tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE_vst.Rdata")
if(!file.exists(tSNEFile)){
  pcsData<-as.data.frame(citeCells@reductions$pca@jackstraw@overall.p.values)
  pcs<-pcsData[pcsData$Score<0.05,"PC"]
  citeCells <- FindNeighbors(object = citeCells, dims = 1:20, verbose = FALSE)
  citeCells <- FindClusters(citeCells, dims.use = 1:20, resolution = 0.6)
  citeCells <- RunTSNE(citeCells, dims.use = 1:20)
  save(citeCells,file=tSNEFile)
} else {load(tSNEFile)}




pcsData<-as.data.frame(citeCells@reductions$pca@jackstraw@overall.p.values)
pcs<-pcsData[pcsData$Score<0.05,"PC"]
citeCells <- RunUMAP(citeCells, dims = pcs)


pdf(paste0(figDir,sampleName,"_",treatment,"_umap.pdf"),width=12,height=8)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", "darkblue","orange",
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F","darkblue","orange")
DimPlot(object = citeCells, reduction = "umap",
         label = TRUE, 
         pt.size = 2, 
         label.size = 6,
         cols = cols
)
dev.off()

citeCells_epithelial<-subset(citeCells,idents=c("0","2","4","15","16"))
pdf(paste0(figDir,sampleName,"_",treatment,"_violin_umap.pdf"),width=12,height=6)
VlnPlot(citeCells_epithelial, ncol=4, features = c("KCTD3", "CLDN3","FOSL2","LAPTM4B","RUSC1","NCKAP1","VTCN1","CALU"))
dev.off()






genes<-c("nFeature_RNA","percent.mt",paste0(c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,sampleName,"_",treatment,"_feature_umap.pdf"),width=12,height=16)
FeaturePlot(object = citeCells, features = genes, ncol = 3, pt.size = 0.5,cols = c("grey", "blue"), 
            reduction = "umap")

dev.off()
#}


set1<-c(
"FARP1"
, "LAMC1"
, "RPL13"
, "KDELR2"
, "KCTD3"
, "ME2"
, "BCL11A"
, "PGRMC1"
, "VEGFA"
)

set2<-c(
"ST3GAL4"
, "CLDN3"
, "ATG4B"
, "KRT19"
, "PTPRF"
, "H2AX"
, "AKT3"
)

set3<-c(
"IGFBP7"
, "HSPG2"
, "DCN"
, "BGN"
, "FOSL2"
, "TIMP2"
, "TIMP3"
, "CXCL9"
, "C1S"
, "C1R"
)

set4<-c(
"RUSC1"
, "KARS"
, "MAP7"
, "ARHGAP30"
, "NCKAP1"
, "RPL13"
, "VTCN1"
, "CALU"
, "PATJ"
, "HACD2"
)

set5<-c(
"LAPTM4B"
, "TXN"
, "PKP3"
, "H1F3"
, "H3FL"
, "RNF44"
)


pdf(paste0(figDir,sampleName,"_",treatment,"_feature_umap_set1.pdf"),width=12,height=9)
FeaturePlot(object = citeCells, features = set1, ncol = 4, pt.size = 0.3,cols = c("grey", "blue"), 
            reduction = "umap")

dev.off()
#}


pdf(paste0(figDir,sampleName,"_",treatment,"_feature_umap_set2.pdf"),width=12,height=6)
FeaturePlot(object = citeCells, features = set2, ncol = 4, pt.size = 0.3,cols = c("grey", "blue"), 
            reduction = "umap")

dev.off()
#}


pdf(paste0(figDir,sampleName,"_",treatment,"_feature_umap_set3.pdf"),width=12,height=9)
FeaturePlot(object = citeCells, features = set3, ncol = 4, pt.size = 0.3,cols = c("grey", "blue"), 
            reduction = "umap")

dev.off()
#}


pdf(paste0(figDir,sampleName,"_",treatment,"_feature_umap_set4.pdf"),width=12,height=9)
FeaturePlot(object = citeCells, features = set4, ncol = 4, pt.size = 0.3,cols = c("grey", "blue"), 
            reduction = "umap")

dev.off()
#}


pdf(paste0(figDir,sampleName,"_",treatment,"_feature_umap_set5.pdf"),width=12,height=3)
FeaturePlot(object = citeCells, features = set5, ncol = 4, pt.size = 0.3,cols = c("grey", "blue"), 
            reduction = "umap")

dev.off()
#}













cite <-cite[,colnames(cite) %in% colnames(citeCells$RNA@data) ]

citeCells[["ADT"]] <- CreateAssayObject(counts = cite)
citeCells <- NormalizeData(citeCells, assay = "ADT", normalization.method = "CLR")
citeCells <- ScaleData(object = citeCells, assay = "ADT")

pdf(paste0(figDir,sampleName,"_",treatment,"_feature_singlet_CITE.pdf"),width=12,height=16)
FeaturePlot(object = citeCells, features = row.names(citeCells$ADT@data), min.cutoff = "q05", max.cutoff = "q90", ncol = 4)
dev.off()

finalObject<-paste0("../project_results/Robjects/",sampleName,"_",treatment,"_CITE.Rdata")
save(citeCells,file=finalObject)

pdf(paste0(figDir,"feature_plot_cite_",sampleName,".pdf"),width=12,height=20)
FeaturePlot(citeCells, features.plot = row.names(cite), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"ridgeplot_",sampleName,".pdf"),width=12,height=4)

RidgePlot(citeCells, features.plot = row.names(cite), cols.use = cols,nCol = 3)

dev.off()
