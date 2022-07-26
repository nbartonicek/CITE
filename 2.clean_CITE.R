library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4497"
sampleName="PBMC_CITE"
treatment="original"
#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))

dataH <- paste0(homeDir,"raw_files/",projectname,"/raw_gene_bc_matrices/GRCh38")
dataM <- paste0(homeDir,"raw_files/",projectname,"/raw_gene_bc_matrices/mm10")

seuratH <- Read10X(dataH)
seuratM <- Read10X(dataM)
seuratX<-rbind(seuratH,seuratM)
seuratMat<-as.matrix(seuratX[,colSums(seuratX)>1000])
#this produces 2646 cells (cutoff 1000)

#species specific
humanIdx<-grepl("GRCh38_",row.names(seuratMat))
speciesRatio<-apply(seuratMat,2,function(x){sum(x[humanIdx])/sum(x)})

#higher proportion of human than mouse
pdf(paste0(figDir,"species_ratio_",sampleName,".pdf"),width=12,height=8)
plot(density(speciesRatio))
dev.off()
pdf(paste0(figDir,"species_ratio_hist_",sampleName,".pdf"),width=12,height=8)
plot(hist(speciesRatio))
dev.off()
#eliminate doublets (0 with cutoff of 0.9)
seuratMatClean<-seuratMat[,speciesRatio>=0.9|speciesRatio<=0.1]

#species list
speciesList<-speciesRatio[speciesRatio>=0.9|speciesRatio<=0.1]
speciesNames<-ifelse(speciesList>0.5,"human","mouse")
speciesColor<-ifelse(speciesList>0.5,"darkblue","red")
save(speciesColor,file="speciesColor.Rdata")
save(speciesRatio,file="speciesRatio.Rdata")

#add the cite data
cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CID4497_PBMC_CITESEQ_DEMUX_S2_hd1_matrix.csv"), 
                     sep = ",", header = TRUE, row.names = 1)
rownames(cite) <- paste0("CITE_", rownames(cite))
cite <- cite[1:(dim(cite)[1]-3),]
cite <-cite[,colnames(cite) %in% colnames(seuratMatClean) ]
cite <-cite[,order(colnames(cite))]
seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(cite) ]


#combine the mouse and human
combined<-CollapseSpeciesExpressionMatrix(seuratMatClean,prefix.1="GRCh38_",prefix.controls="mm10_",features.controls.toKeep=500)

citeCells <- CreateSeuratObject(
  raw.data = combined, 
  min.cells = 3, 
  min.genes = 100, 
  project = "4497_cite"
)

save(citeCells,file=paste0("../project_results/Robjects/",sampleName,"_earlySeurat.Rdata"))

mito.genes <- grep(pattern = "^MT-|^mm10___mt-", x = rownames(x = citeCells@data), value = TRUE)
percent.mito <- Matrix::colSums(citeCells@raw.data[mito.genes, ])/Matrix::colSums(citeCells@raw.data)

citeCells <- AddMetaData(object = citeCells, metadata = percent.mito, col.name = "percent.mito")
citeCells <- FilterCells(object = citeCells, subset.names = c("nGene", "percent.mito"), low.thresholds = c(5, -Inf), high.thresholds = c(7000, 0.3))

#citeCells <- AddMetaData(object = citeCells, metadata = speciesNames, col.name = "species")
pdf(paste0(figDir,"vlnplot_",sampleName,".pdf"),width=12,height=8)
VlnPlot(object = citeCells, 
        features.plot = c("nGene", "nUMI", "percent.mito"), 
        nCol = 3)
dev.off()
pdf(paste0(figDir,"stats_",sampleName,".pdf"),width=12,height=4)

par(mfrow = c(1, 3))

plot(citeCells@meta.data$nUMI,citeCells@meta.data$nGene,pch=15,cex=0.6,col=speciesColor)
plot(citeCells@meta.data$nUMI,citeCells@meta.data$percent.mito,pch=15,cex=0.6,col=speciesColor)
plot(citeCells@meta.data$nGene,citeCells@meta.data$percent.mito,pch=15,cex=0.6,col=speciesColor)
dev.off()

citeCells <- NormalizeData(object = citeCells)

citeCells <- FindVariableGenes(citeCells, do.plot = T, y.cutoff = 0.5)

citeCells <- ScaleData(citeCells, display.progress = FALSE)

citeCells <- RunPCA(citeCells, pcs.print = 0)
pdf(paste0(figDir,"PCA_bowl_plot_",sampleName,".pdf"),width=12,height=8)
PCElbowPlot(citeCells)
dev.off()
citeCells <- FindClusters(citeCells, dims.use = 1:15, print.output = FALSE)
citeCells <- RunTSNE(citeCells, dims.use = 1:15)
save(citeCells,file=paste0("../project_results/Robjects/",sampleName,"_midSeurat.Rdata"))

pdf(paste0(figDir,"TSNE_",sampleName,".pdf"),width=12,height=8)

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
temp <- FindAllMarkers(citeCells, max.cells.per.ident = 100, logfc.threshold = log(2), 
                                   only.pos = TRUE, min.diff.pct = 0.3, do.print = F)

cluster2.markers <- FindMarkers(object = citeCells, ident.1 = 2, ident.2 = 1,thresh.use = 0.25, 
    test.use = "roc", only.pos = TRUE)

predef_genes_human<-c("VCL")

#pdf(paste0(figDir,sampleName,"_tSNE_with_predefined_genes_human.pdf"),width=8,height=8)
pdf(paste0(figDir,"featurePlot_",sampleName,".pdf"),width=12,height=8)
FeaturePlot(object = citeCells, features.plot = predef_genes_human, cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
dev.off()
cite <-cite[colnames(cite) %in% colnames(citeCells@data) ]

citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)

pdf(paste0(figDir,"feature_plot_cite_",sampleName,".pdf"),width=12,height=20)
FeaturePlot(citeCells, features.plot = row.names(cite), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"cross_cite_",sampleName,".pdf"),width=12,height=4)

par(mfrow = c(1, 3))
GenePlot(citeCells, gene1 = "CITE_H1", gene2 = "CITE_H2", cex.use = 0.5, col.use = cols)
GenePlot(citeCells, gene1 = "CITE_H1", gene2 = "CITE_H3", cex.use = 0.5, col.use = cols)
GenePlot(citeCells, gene1 = "CITE_H2", gene2 = "CITE_H3", cex.use = 0.5, col.use = cols)

dev.off()

pdf(paste0(figDir,"ridgeplot_",sampleName,".pdf"),width=12,height=4)

RidgePlot(citeCells, features.plot = c("CITE_H1", "CITE_H2", "CITE_H3"), cols.use = cols,nCol = 3)

dev.off()



#testing for dead cells
FeaturePlot(citeCells, features.plot = c("MT-ND1","mm10___mt-Nd1","MYC","MGAT1"), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 3, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
save(citeCells,file="../project_results/Robjects/4497_PBMC_CITE.Rdata")


