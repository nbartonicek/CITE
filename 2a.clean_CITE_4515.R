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
cite <-cite[,order(colnames(cite))]


earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_earlySeurat.Rdata")
if(!file.exists(earlySeuratFile)){
  dataH <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/GRCh38")
  dataM <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/mm10")

  seuratH <- Read10X(dataH)
  #seuratM <- Read10X(dataM)
  #seuratX<-rbind(seuratH,seuratM)
  seuratX <- seuratH
  seuratMat<-as.matrix(seuratX[,colSums(seuratX)>1000])
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
    raw.data = seuratMatClean, 
    min.cells = 3, 
    min.genes = 100, 
    project = sampleName
  )

  save(citeCells,file=earlySeuratFile)
} else {load(earlySeuratFile)}

mito.genes <- grep(pattern = "^MT-|^mm10___mt-", x = rownames(x = citeCells@data), value = TRUE)
percent.mito <- Matrix::colSums(citeCells@raw.data[mito.genes, ])/Matrix::colSums(citeCells@raw.data)

citeCells <- AddMetaData(object = citeCells, metadata = percent.mito, col.name = "percent.mito")
citeCells <- FilterCells(object = citeCells, subset.names = c("nGene", "percent.mito"), low.thresholds = c(250, -Inf), high.thresholds = c(Inf, 0.1))

#citeCells <- AddMetaData(object = citeCells, metadata = speciesNames, col.name = "species")
pdf(paste0(figDir,"vlnplot_",sampleName,".pdf"),width=12,height=8)
VlnPlot(object = citeCells, 
        features.plot = c("nGene", "nUMI", "percent.mito"), 
        nCol = 3)
dev.off()
pdf(paste0(figDir,"stats_",sampleName,".pdf"),width=12,height=4)

par(mfrow = c(1, 3))

plot(citeCells@meta.data$nUMI,citeCells@meta.data$nGene,pch=15,cex=0.6)
plot(citeCells@meta.data$nUMI,citeCells@meta.data$percent.mito,pch=15,cex=0.6)
plot(citeCells@meta.data$nGene,citeCells@meta.data$percent.mito,pch=15,cex=0.6)
dev.off()

citeCells <- NormalizeData(object = citeCells)
citeCells <- FindVariableGenes(citeCells, do.plot = T, y.cutoff = 0.5)
citeCells <- ScaleData(citeCells, display.progress = FALSE)
citeCells <- RunPCA(citeCells, pcs.print = 0,pcs.compute = 50)
citeCells <- ProjectPCA(object = citeCells, pcs.store = 50, do.print = FALSE)

pdf(paste0(figDir,"PCA_heatmap_plot_",sampleName,".pdf"),width=16,height=30)
PCHeatmap(object = citeCells, pc.use = 1:30, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
dev.off()

citeCells <- JackStraw(object = citeCells, num.pc = 50,num.replicate = 100, display.progress = FALSE)
pdf(paste0(figDir,"PCA_jackstraw_plot_",sampleName,".pdf"),width=12,height=40)
JackStrawPlot(object = citeCells, PCs = 1:100)
dev.off()

#save(citeCells,file=paste0("../project_results/Robjects/",sampleName,"_jackstraw.Rdata"))
#break()

pdf(paste0(figDir,"PCA_bowl_plot_",sampleName,".pdf"),width=12,height=8)
PCElbowPlot(citeCells)
dev.off()

citeCells <- FindClusters(citeCells, dims.use = 1:25, resolution = 0.6,print.output = FALSE)
citeCells <- RunTSNE(citeCells, dims.use = 1:25,perplexity=10)

FindMarkers(object = citeCells, ident.1 = 1, min.pct = 0.25)


save(citeCells,file=paste0("../project_results/Robjects/",sampleName,"_midSeurat.Rdata"))


pdf(paste0(figDir,"TSNE_",sampleName,"4_25_p10.pdf"),width=12,height=8)

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

genes<-c("nGene","percent.mito","CD14","EPCAM","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8", "PECAM1", "CD34", "ERBB2")


#pdf(paste0(figDir,sampleName,"_tSNE_with_predefined_genes_human.pdf"),width=8,height=8)
pdf(paste0(figDir,"featurePlot_",sampleName,"4_25_p10.pdf"),width=12,height=8)
FeaturePlot(object = citeCells, features.plot = genes, cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
dev.off()



#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)

pdf(paste0(figDir,"feature_plot_cite_",sampleName,".pdf"),width=12,height=20)
FeaturePlot(citeCells, features.plot = row.names(cite), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"ridgeplot_",sampleName,".pdf"),width=12,height=4)

RidgePlot(citeCells, features.plot = row.names(cite), cols.use = cols,nCol = 3)

dev.off()


genes<-c("CD14","CD16","COL1A1","EPCAM","CD19","CD3","","","")
#testing for dead cells
FeaturePlot(citeCells, features.plot = c("MT-ND1","mm10___mt-Nd1","MYC","MGAT1"), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 3, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
save(citeCells,file="../project_results/Robjects/4497_PBMC_CITE.Rdata")

load(paste0("../project_results/Robjects/4515_CITE_annotated_imputed.Rdata"))

pdf(paste0(figDir,"TSNE_",sampleName,"annotated_imputed.pdf"),width=12,height=8)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F","pink","black")
TSNEPlot(object = annotated_citeCells, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         colors.use = cols
)
dev.off()


