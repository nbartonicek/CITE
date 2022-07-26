library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)

#Define project and sample/methodology within it within it
projectname="3921_3963_4066"
sampleName="3921_3963_4066_HTO"
treatment="v3"

#directory structure
homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))


#load CITE or HASH data
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


#create Seurat object and save intermediary steps for faster re-runs 
earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_",treatment,"_earlySeurat.Rdata")

if(!file.exists(earlySeuratFile)){
  dataH <- paste0(homeDir,"raw_files/",projectname,"/outs/filtered_gene_bc_matrices/GRCh38")
  seuratX <- Read10X(dataH)
  
  #species specific - calculate ratio of species if multiple
  #humanIdx<-grepl("GRCh38_",row.names(seuratX))
  #speciesRatio<-apply(seuratX,2,function(x){sum(x[humanIdx])/sum(x)})
  #seuratMatClean<-seuratX[,speciesRatio>=0.9|speciesRatio<=0.1]


  #if only one species
  seuratMatClean<-seuratX
  
  #take only those cells that have cite and/or hash data if there
  #seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(hash) ]
  #seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(cite) ]


  #load the Seurat object, a gene should be present in at least 3 cells
  #and only cells with >100 genes should be accepted
  citeCells <- CreateSeuratObject(
    counts = seuratMatClean, 
    min.cells = 3, 
    min.features = 100, 
    project = sampleName
  )

  save(citeCells,file=earlySeuratFile)
} else {load(earlySeuratFile)}

#calculate percentage of mitochondrial genes
citeCells <- PercentageFeatureSet(object = citeCells, pattern = "^MT-", col.name = "percent.mt")

#plot violin plot for QC
pdf(paste0(figDir,"vlnplotS3_",sampleName,".pdf"),width=12,height=8)
VlnPlot(object = citeCells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#eliminate cells with high MITO content
citeCells <- subset(x = citeCells, percent.mt < 25)

#normalize and find principal components
sctranscormFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_sctransform.Rdata")
if(!file.exists(sctranscormFile)){
  citeCells <- SCTransform(object = citeCells, verbose = FALSE)
  citeCells <- RunPCA(object = citeCells, verbose = FALSE,npcs = 100)
  citeCells <- JackStraw(object = citeCells, dims=100,num.replicate=100,verbose=T)
  citeCells <- ScoreJackStraw(citeCells,dims=1:100)
  save(citeCells,file=sctranscormFile)
} else {load(sctranscormFile)}

#make tsne
tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE_vst.Rdata")
if(!file.exists(tSNEFile)){
  pcsData<-as.data.frame(citeCells@reductions$pca@jackstraw@overall.p.values)
  pcs<-pcsData[pcsData$Score<0.05,"PC"]
  citeCells <- FindNeighbors(object = citeCells, dims = pcs, verbose = FALSE)
  citeCells <- FindClusters(citeCells, dims.use = pcs, resolution = 0.6)
  citeCells <- RunTSNE(citeCells, dims.use = pcs)
  save(citeCells,file=tSNEFile)
} else {load(tSNEFile)}


#plot tsne
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


#plot classic genes
genes<-c("nFeature_RNA","percent.mt",paste0(c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,sampleName,"_",treatment,"_feature_vst.pdf"),width=12,height=16)
FeaturePlot(object = citeCells, features = genes, ncol = 3, pt.size = 0.5,cols = c("grey", "blue"), 
            reduction = "tsne")

dev.off()


#add hash data
#hash <-hash[,colnames(hash) %in% colnames(citeCells$RNA@data) ]
#citeCells[["HTO"]] <- CreateAssayObject(counts = hash)
#citeCells <- NormalizeData(object = citeCells, assay = "HTO", normalization.method = "CLR")
#citeCells <- HTODemux(citeCells, assay = "HTO", positive.quantile = 0.99)
#table(citeCells$HTO_classification.global)
#
#Idents(object = citeCells) <- "HTO_maxID"
#pdf(paste0(figDir,"ridgeplot_hashv3_",sampleName,".pdf"),width=12,height=4)
#RidgePlot(object = citeCells, assay = "HTO", features = rownames(x = citeCells[["HTO"]])[1:3], 
#    ncol = 3)
#dev.off()
#
#Idents(object = citeCells) <- "HTO_classification.global"
#citeCells.singlet <- subset(x = citeCells, idents = "Singlet")
#
#add cite data
#cite <-cite[,colnames(cite) %in% colnames(citeCells.singlet$RNA@data) ]
#citeCells.singlet[["ADT"]] <- CreateAssayObject(counts = cite)
#citeCells.singlet <- NormalizeData(citeCells.singlet, assay = "ADT", normalization.method = "CLR")
#citeCells.singlet <- ScaleData(object = citeCells.singlet, assay = "ADT")
#
#pdf(paste0(figDir,sampleName,"_",treatment,"_feature_singlet_CITE.pdf"),width=12,height=16)
#FeaturePlot(object = citeCells.singlet, features = row.names(citeCells.singlet$ADT@data), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)
#dev.off()

