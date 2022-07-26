library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
#library(DropletUtils)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="3921_3963_4066"
sampleName="3921_3963_4066_HTO"
treatment="v2"

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

earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_",treatment,"_earlySeurat.Rdata")
if(!file.exists(earlySeuratFile)){
  dataH <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/GRCh38")

  seuratX <- Read10X(dataH)
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


  #combine the mouse and human
  
  citeCells <- CreateSeuratObject(
    raw.data = seuratMatClean, 
    min.cells = 3, 
    min.genes = 100, 
    project = sampleName
  )

  save(citeCells,file=earlySeuratFile)
} else {load(earlySeuratFile)}


jackstrawFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_jackstraw.Rdata")
if(!file.exists(jackstrawFile)){

  mito.genes <- grep(pattern = "^MT-|^mm10___mt-", x = rownames(x = citeCells@data), value = TRUE)
  percent.mito <- Matrix::colSums(citeCells@raw.data[mito.genes, ])/Matrix::colSums(citeCells@raw.data)

  citeCells <- AddMetaData(object = citeCells, metadata = percent.mito, col.name = "percent.mito")
  citeCells <- FilterCells(object = citeCells, subset.names = c("nGene", "percent.mito","nUMI"), low.thresholds = c(100, -Inf, 500), high.thresholds = c(Inf, 0.25,Inf))
  citeCells <- NormalizeData(object = citeCells)
  citeCells <- FindVariableGenes(citeCells, do.plot = T,  mean.function = ExpMean,dispersion.function = LogVMR, x.low.cutoff = 0.125, x.high.cutoff = 3,y.cutoff = 0.5)
  citeCells <- ScaleData(citeCells, vars.to.regress=c("nUMI","nGene"),display.progress = FALSE)
  citeCells <- RunPCA(citeCells, pcs.print = 0,npcs = 100)
  citeCells <- ProjectPCA(object = citeCells, pcs.store = 100, do.print = FALSE)
  citeCells <- JackStraw(object = citeCells, num.pc = 100,num.replicate = 100)
  citeCells <- JackStrawPlot(object = citeCells, PCs = 1:20)

  save(citeCells,file=jackstrawFile)
}else{load(jackstrawFile)}
pdf(paste0(figDir,"vlnplot_",sampleName,".pdf"),width=12,height=6)
VlnPlot(object = citeCells, features = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE.Rdata")
if(!file.exists(tSNEFile)){


  pcsData<-as.data.frame(citeCells@dr$pca@jackstraw@overall.p.values)
  pcs<-pcsData[pcsData$Score<0.01,"PC"]
  #for(perplexity in perplexities){
  citeCells <- FindClusters(citeCells, dims.use = pcs, resolution = resolution,print.output = FALSE)
  citeCells <- RunTSNE(citeCells, dims.use = pcs, perplexity=perplexity)
  save(citeCells,file=tSNEFile)
} else {load(tSNEFile)}

dev.off()


#citeCells=UpdateSeuratObject(citeCells)

pdf(paste0(figDir,sampleName,"_",treatment,"_tsne.pdf"),width=12,height=8)

cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F")
TSNEPlot(object = citeCells, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         colors.use = cols
)

dev.off()

#genes<-c("nFeature_RNA","percent.mito","mm10---Rps15a","ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")
#genes<-c("nGene","percent.mito","ACTA2","CD34","PECAM1","PTPRC","CD68","FCGR3A","CD19","CD3D","FOXP3")

genes<-c("nGene","percent.mito",paste0(c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,sampleName,"_",treatment,"_feature.pdf"),width=12,height=16)
FeaturePlot(object = citeCells, features.plot = genes, nCol = 4, pt.size = 0.5,cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

dev.off()
#}


#pdf(paste0(figDir,sampleName,"_tSNE_with_predefined_genes_human.pdf"),width=8,height=8)
pdf(paste0(figDir,"featurePlot_",sampleName,"_",treatment,".pdf"),width=12,height=8)
FeaturePlot(object = citeCells, pt.size = 0.2,    features = genes, cols = c("grey", "blue"), 
            reduction = "tsne",label.size = 0.5)
dev.off()


#add CITE data
hash <-hash[,colnames(hash) %in% colnames(citeCells@data) ]


citeCells <- SetAssayData(citeCells, assay.type = "HTO", slot = "raw.data", new.data = hash)
citeCells <- NormalizeData(citeCells, assay.type = "HTO", normalization.method = "genesCLR")
citeCells <- HTODemux(citeCells,assay.type = "HTO",positive_quantile =  0.99, print.output = TRUE)

citeCells <- SetAllIdent(citeCells,"hto_classification_global")

pdf(paste0(figDir,"feature_plot_hash2_max3_",sampleName,".pdf"),width=12,height=4)
FeaturePlot(citeCells, features.plot = c("HTO1", "HTO2", "HTO3"), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 3, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()

pdf(paste0(figDir,"ridgeplot_",sampleName,".pdf"),width=12,height=4)
RidgePlot(citeCells,features.plot = rownames(GetAssayData(citeCells,assay.type = "HTO"))[1:3],nCol = 3)
dev.off()

citeCells <- SetAllIdent(citeCells,"hto_classification")
cols<-c("#b2182b", "pink", "purple", "#053061", 
                         "orange", "#1c9099", 
                        "gray")

pdf(paste0(figDir,"tsneplot_",sampleName,".pdf"),width=12,height=8)
TSNEPlot(citeCells,group.by = "hto_classification", 
         pt.size = 1, 
         label.size = 6,
         colors.use = cols)
dev.off()
citeCells <- SetAllIdent(citeCells,"hto_classification_global")

pdf(paste0(figDir,"hto_classification_global_",sampleName,".pdf"),width=12,height=8)
TSNEPlot(citeCells,group.by = "hto_classification", 
         pt.size = 1, 
         label.size = 6,
         colors.use = cols)
dev.off()
