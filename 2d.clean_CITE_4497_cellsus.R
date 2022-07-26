library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
#library(DropletUtils)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4497-2_CELLSUS"
sampleName="4497-2_CELLSUS_CITE"
treatment="original"
#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))

cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,".tsv"), 
                       sep = ",", header = TRUE, row.names = 1)

cite <- cite[!(row.names(cite) %in% c("bad_struct","no_match","total_reads")),]

rownames(cite) <- paste0("CITE_", rownames(cite))
cite <-cite[,order(colnames(cite))]


earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_earlySeurat.Rdata")
if(!file.exists(earlySeuratFile)){
  dataH <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/GRCh38")
  dataM <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/mm10")

  seuratH <- Read10X(dataH)
  seuratM <- Read10X(dataM)
  seuratX<-rbind(seuratH,seuratM)
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
  seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(cite) ]


  #combine the mouse and human
  combined<-CollapseSpeciesExpressionMatrix(seuratMatClean,prefix="GRCh38_",controls="mm10_",ncontrols=500)

  citeCells <- CreateSeuratObject(
    counts = combined, 
    min.cells = 3, 
    min.features = 100, 
    project = sampleName
  )

  save(citeCells,file=earlySeuratFile)
} else {load(earlySeuratFile)}


mito.genes <- grep(pattern = "^MT-|^mm10---mt-", x = rownames(x = citeCells@assays$RNA@data), value = TRUE)
percent.mito <- Matrix::colSums(citeCells@assays$RNA@counts[mito.genes, ])/Matrix::colSums(citeCells@assays$RNA@counts)
citeCells <- AddMetaData(object = citeCells, metadata = percent.mito, col.name = "percent.mito")

pdf(paste0(figDir,"vlnplotS3_",sampleName,".pdf"),width=12,height=6)
VlnPlot(object = citeCells, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

citeCells <- subset(x = citeCells, subset = percent.mito < 0.25)
#from 3124 to 4772 with Seurat3

pdf(paste0(figDir,"statsS3_",sampleName,".pdf"),width=12,height=4)
par(mfrow = c(1, 3))
plot(citeCells@meta.data$nCount_RNA,citeCells@meta.data$nFeature_RNA,pch=15,cex=0.6)
plot(citeCells@meta.data$nCount_RNA,citeCells@meta.data$percent.mito,pch=15,cex=0.6)
plot(citeCells@meta.data$nFeature_RNA,citeCells@meta.data$percent.mito,pch=15,cex=0.6)
dev.off()

citeCells <- NormalizeData(object = citeCells, normalization.method = "LogNormalize", scale.factor = 10000)
citeCells <- FindVariableFeatures(citeCells, do.plot = T, selection.method = "mean.var.plot", mean.cutoff = c(0.1, 8), dispersion.cutoff = c(0.5, Inf))


pdf(paste0(figDir,"variableGene_",sampleName,".pdf"),width=8,height=4)
VariableFeaturePlot(object = citeCells)
dev.off()


#length(x = VariableFeatures(object = citeCells))
#3995 features... seems a lot

jackstrawFile=paste0("../project_results/Robjects/",sampleName,"_jackstrawS3.Rdata")

citeCells <- ScaleData(citeCells, vars.to.regress=c("nCount_RNA","nFeature_RNA"),display.progress = FALSE)
citeCells <- RunPCA(citeCells, pcs.print = 0,npcs = 100)
citeCells <- ProjectDim(object = citeCells)
citeCells <- JackStraw(object = citeCells, dims = 50,num.replicate = 100)
citeCells <- ScoreJackStraw(object = citeCells, dims = 1:50)

save(citeCells,file=jackstrawFile)


pdf(paste0(figDir,"PCA_jackstraw_plot_",sampleName,".pdf"),width=12,height=40)
JackStrawPlot(object = citeCells, dims = 1:50)
dev.off()

pdf(paste0(figDir,"PCA_heatmap_plot_",sampleName,".pdf"),width=16,height=30)
PCHeatmap(object = citeCells, pc.use = 1:30, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
dev.off()

pdf(paste0(figDir,"PCA_bowl_plot_",sampleName,".pdf"),width=12,height=8)
ElbowPlot(citeCells)
dev.off()

pcsData<-as.data.frame(citeCells@reductions$pca@jackstraw@overall.p.values)
pcs<-pcsData[pcsData$Score<0.01,"PC"]
  #for(perplexity in perplexities){

tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_tSNE.Rdata")
if(!file.exists(tSNEFile)){
  citeCells <- FindNeighbors(citeCells, dims=pcs)
  citeCells <- FindClusters(citeCells, resolution = 0.6)
  citeCells <- RunTSNE(citeCells, dims = pcs)
  save(citeCells,file=tSNEFile)
} else {load(tSNEFile)}  



pdf(paste0(figDir,"TSNE_",sampleName,".pdf"),width=12,height=8)

cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F")
DimPlot(object = citeCells, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         cols = cols
)
dev.off()

genes<-c("nFeature_RNA","percent.mito","mm10---Rps15a","ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")
#genes<-c("nGene","percent.mito","ACTA2","CD34","PECAM1","PTPRC","CD68","FCGR3A","CD19","CD3D","FOXP3")


#pdf(paste0(figDir,sampleName,"_tSNE_with_predefined_genes_human.pdf"),width=8,height=8)
pdf(paste0(figDir,"featurePlot_",sampleName,".pdf"),width=12,height=8)
FeaturePlot(object = citeCells, pt.size = 0.2,    features = genes, cols = c("grey", "blue"), 
            reduction = "tsne",label.size = 0.5)
dev.off()



#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells@assays$RNA@data) ]
row.names(cite)=gsub(" ","",row.names(cite))
citeCells <- SetAssayData(citeCells, assay.type = "ADT", slot = "counts", new.data = as.matrix(cite))
citeCells <- NormalizeData(citeCells, assay.type = "ADT", normalization.method = "CLR")
#citeCells <- ScaleData(citeCells, assay.type = "ADT", display.progress = FALSE)

pdf(paste0(figDir,"feature_plot_cite_",sampleName,"1.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[1:20], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"2.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[21:40], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"3.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[41:60], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"4.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[61:80], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"5.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[81:93], min.cutoff = "q05", max.cutoff = "q95", 
    ncol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"ridgeplot_",sampleName,".pdf"),width=12,height=4)

RidgePlot(citeCells, features = row.names(cite), cols = cols,ncol = 3)

dev.off()




###### do the clustering based on CITE #####

DefaultAssay(object = citeCells) <- "CITE"

citeCells_cite <- RunPCA(citeCells, pc.genes = rownames(cite), assay.type = "CITE", 
    pcs.print = 0)


pdf(paste0(figDir,"PCAplot_",sampleName,".pdf"),width=6,height=4)

DimPlot(citeCells_cite, pt.size = 0.5)

dev.off()












