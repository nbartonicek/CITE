library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
#library(DropletUtils)

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



cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_CITE.tsv"), 
                       sep = ",", header = TRUE, row.names = 1)

cite <- cite[!(row.names(cite) %in% c("bad_struct","no_match","total_reads")),]

#rownames(hash) <- paste0("HTO_", rownames(hash))
cite <-cite[,order(colnames(cite))]



earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_",treatment,"_earlySeurat.Rdata")
if(!file.exists(earlySeuratFile)){
  dataH <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/GRCh38")

  seuratX <- Read10X(dataH)
  #seuratMat<-as.matrix(seuratX[,colSums(seuratX)>500])
  #this produces 4153 cells (cutoff 500)
  load("noncancer.Rdata")
  #whitelist=read.table(file=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_whitelist.txt"))
  whitelist=noncancer
  seuratMat=seuratX[,colnames(seuratX) %in% whitelist]
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




tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE_vst.Rdata")
if(!file.exists(tSNEFile)){
  pcsData<-as.data.frame(citeCells@reductions$pca@jackstraw@overall.p.values)
  pcs<-pcsData[pcsData$Score<0.05,"PC"]
  citeCells <- FindNeighbors(object = citeCells, dims = pcs, verbose = FALSE)
  citeCells <- FindClusters(citeCells, dims.use = pcs, resolution = 0.6)
  citeCells <- RunTSNE(citeCells, dims.use = pcs)
  save(citeCells,file=tSNEFile)
} else {load(tSNEFile)}

dev.off()


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
citeCells <- HTODemux(citeCells, assay = "HTO", positive.quantile = 0.99)
table(citeCells$HTO_classification.global)

pdf(paste0(figDir,"umap_HTO_classification_",sampleName,"_",treatment,".pdf"),width=10,height=8)
DimPlot(object = citeCells.singlet, group.by = "HTO_classification",pt.size = 2,cols=cols)
dev.off()

Idents(object = citeCells) <- "HTO_maxID"
pdf(paste0(figDir,"ridgeplot_hash_",sampleName,"_",treatment,".pdf"),width=12,height=4)
RidgePlot(object = citeCells, assay = "HTO", features = rownames(x = citeCells[["HTO"]])[1:3], 
    ncol = 3)
dev.off()

Idents(object = citeCells) <- "HTO_classification.global"
citeCells.singlet <- subset(x = citeCells, idents = "Singlet")

# Calculate tSNE embeddings with a distance matrix
citeCells.singlet <- FindVariableFeatures(object = citeCells.singlet, selection.method = "mean.var.plot")
citeCells.singlet <- ScaleData(object = citeCells.singlet, features = VariableFeatures(object = citeCells.singlet))
citeCells.singlet <- RunPCA(object = citeCells.singlet, features = VariableFeatures(object = citeCells.singlet))
citeCells.singlet <- FindNeighbors(object = citeCells.singlet, reduction = "pca", dims = 1:20)
citeCells.singlet <- FindClusters(object = citeCells.singlet, resolution = 0.6)
citeCells.singlet <- RunTSNE(object = citeCells.singlet, reduction = "pca", dims = 1:20)

cols<-c("#b2182b", "#053061", "#1c9099")
pdf(paste0(figDir,"tsne_singlets_hashv3_",sampleName,".pdf"),width=10,height=8)
DimPlot(object = citeCells.singlet, group.by = "HTO_classification",pt.size = 2,cols=cols)
dev.off()


pdf(paste0(figDir,sampleName,"_",treatment,"_tsne_singlet.pdf"),width=12,height=8)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F","darkblue","orange")
DimPlot(object = citeCells.singlet, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         cols = cols
)
dev.off()


genes<-c("nFeature_RNA","percent.mt",paste0(c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,sampleName,"_",treatment,"_feature_singlet.pdf"),width=12,height=16)
FeaturePlot(object = citeCells.singlet, features = genes, ncol = 3, pt.size = 0.5,cols = c("grey", "blue"), 
            reduction = "tsne")

dev.off()
#}

cite <-cite[,colnames(cite) %in% colnames(citeCells.singlet$RNA@data) ]

citeCells.singlet[["ADT"]] <- CreateAssayObject(counts = cite)
citeCells.singlet <- NormalizeData(citeCells.singlet, assay = "ADT", normalization.method = "CLR")
citeCells.singlet <- ScaleData(object = citeCells.singlet, assay = "ADT")

pdf(paste0(figDir,sampleName,"_",treatment,"_feature_singlet_CITE.pdf"),width=12,height=16)
FeaturePlot(object = citeCells.singlet, features = row.names(citeCells.singlet$ADT@data), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)
dev.off()

