library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4497"
sample="CID4497-PDX_Hash"
treatment="original"
treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sample,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))

dataH <- paste0(homeDir,"raw_files/",projectname,"/",sample,"/raw_gene_bc_matrices/GRCh38")
dataM <- paste0(homeDir,"raw_files/",projectname,"/",sample,"/raw_gene_bc_matrices/mm10")

seuratH <- Read10X(dataH)
seuratM <- Read10X(dataM)
seuratX<-rbind(seuratH,seuratM)
seuratMat<-as.matrix(seuratX[,colSums(seuratX)>1000])
#this produces 5330 cells (cutoff 1000)

#species specific
humanIdx<-grepl("GRCh38_",row.names(seuratMat))
speciesRatio<-apply(seuratMat,2,function(x){sum(x[humanIdx])/sum(x)})

#higher proportion of human than mouse
pdf(paste0(figDir,"species_ratio_",sample,".pdf"),width=12,height=8)
plot(density(speciesRatio))
dev.off()
pdf(paste0(figDir,"species_ratio_hist_",sample,".pdf"),width=12,height=8)
plot(hist(speciesRatio))
dev.off()
#eliminate doublets (only 6 with cutoff of 0.9)
seuratMatClean<-seuratMat[,speciesRatio>=0.9|speciesRatio<=0.1]

#species list
speciesList<-speciesRatio[speciesRatio>=0.9|speciesRatio<=0.1]
speciesNames<-ifelse(speciesList>0.5,"human","mouse")
speciesColor<-ifelse(speciesList>0.5,"darkblue","red")

#add the HASH data
hash <- read.csv(paste0(homeDir,"raw_files/",projectname,"/Human_HASHSEQ_S1_hd2_F.csv"), 
                     sep = ",", header = TRUE, row.names = 1)
rownames(hash) <- paste0("HTO_", rownames(hash))
hash <- hash[1:(dim(hash)[1]-3),]
hash <-hash[colnames(hash) %in% colnames(seuratMatClean) ]
hash <-hash[,order(colnames(hash))]
seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(hash) ]

save(hash,file="../project_results/Robjects/Human_HASHSEQ_S1_hd2_F_clean.Rdata")

#combine the mouse and human
combined<-CollapseSpeciesExpressionMatrix(seuratMatClean,prefix.1="GRCh38_",prefix.controls="mm10_",features.controls.toKeep=500)

hashedCells <- CreateSeuratObject(
  raw.data = combined, 
  min.cells = 3, 
  min.genes = 100, 
  project = "4497_hash"
)


mito.genes <- grep(pattern = "^MT-|^mm10___mt-", x = rownames(x = hashedCells@data), value = TRUE)
percent.mito <- Matrix::colSums(hashedCells@raw.data[mito.genes, ])/Matrix::colSums(hashedCells@raw.data)

hashedCells <- AddMetaData(object = hashedCells, metadata = percent.mito, col.name = "percent.mito")
hashedCells <- FilterCells(object = hashedCells, subset.names = c("nGene", "nUMI","percent.mito"), low.thresholds = c(200, 500,-Inf), high.thresholds = c(Inf, Inf,0.3))

#hashedCells <- AddMetaData(object = hashedCells, metadata = speciesNames, col.name = "species")
pdf(paste0(figDir,"vlnplot_",sample,".pdf"),width=12,height=8)
VlnPlot(object = hashedCells, 
        features.plot = c("nGene", "nUMI", "percent.mito"), 
        nCol = 3)
dev.off()
pdf(paste0(figDir,"stats_",sample,".pdf"),width=12,height=4)

par(mfrow = c(1, 3))
#GenePlot(object = hashedCells, 
#         gene1 = "nUMI", 
#         gene2 = "percent.mito",
#         cex.use = 0.4)
#GenePlot(object = hashedCells, 
#         gene1 = "nUMI", 
#         gene2 = "nGene",
#         cex.use = 0.4)
#GenePlot(object = hashedCells, 
#         gene1 = "nGene", 
#         gene2 = "percent.mito", 
#         cex.use = 0.4)
plot(hashedCells@meta.data$nUMI,hashedCells@meta.data$nGene,pch=15,cex=0.6,col=speciesColor)
plot(hashedCells@meta.data$nUMI,hashedCells@meta.data$percent.mito,pch=15,cex=0.6,col=speciesColor)
plot(hashedCells@meta.data$nGene,hashedCells@meta.data$percent.mito,pch=15,cex=0.6,col=speciesColor)
dev.off()

hashedCells <- NormalizeData(object = hashedCells)

hashedCells <- FindVariableGenes(hashedCells, do.plot = T, y.cutoff = 0.5)

hashedCells <- ScaleData(hashedCells, vars.to.regress=c("nUMI","nGene"),display.progress = FALSE)

hashedCells <- RunPCA(hashedCells, pcs.print = 0)
pdf(paste0(figDir,"PCA_bowl_plot_",sample,".pdf"),width=12,height=8)
PCElbowPlot(hashedCells)
dev.off()
hashedCells <- FindClusters(hashedCells, dims.use = 1:15, print.output = FALSE)
hashedCells <- RunTSNE(hashedCells, dims.use = 1:15)
pdf(paste0(figDir,"TSNE_",sample,".pdf"),width=12,height=8)

cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F")
TSNEPlot(object = hashedCells, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         colors.use = cols
)
dev.off()


#do the hashing
joint_bcs <- intersect(colnames(hashedCells@data),colnames(hash))
#hashedCells_sparse <- hashedCells@data[,joint_bcs]
hash <- as.matrix(hash[,joint_bcs])

hashedCells <- SetAssayData(hashedCells, assay.type = "HTO", slot = "raw.data", new.data = hash)
hashedCells <- NormalizeData(hashedCells, assay.type = "HTO", normalization.method = "genesCLR")
#hashedCells <- ScaleData(hashedCells, assay.type = "HTO", display.progress = FALSE)
hashedCells <- HTODemux(hashedCells,assay.type = "HTO",positive_quantile =  0.99,print.output = FALSE)
print (table(hashedCells@meta.data$hto_classification_global))

hashedCells <- SetAllIdent(hashedCells,id = "hash_maxID")
pdf(paste0(figDir,"ridgeplot2_",sample,".pdf"),width=12,height=4)
RidgePlot(hashedCells,features.plot = rownames(GetAssayData(hashedCells,assay.type = "HTO"))[1:3],nCol = 3)
dev.off()

pdf(paste0(figDir,"feature_plot_hash2_",sample,".pdf"),width=12,height=4)
FeaturePlot(hashedCells, features.plot = c("HTO_H1", "HTO_H2", "HTO_H3"), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 3, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_hash2_max_",sample,".pdf"),width=12,height=4)
FeaturePlot(hashedCells, features.plot = c("HTO_H1", "HTO_H2", "HTO_H3"), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 3, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"cross_hash2_",sample,".pdf"),width=8,height=4)

par(mfrow = c(1, 3))
GenePlot(hashedCells, gene1 = "HTO_H1", gene2 = "HTO_H2", cex.use = 0.5, col.use = cols)
GenePlot(hashedCells, gene1 = "HTO_H1", gene2 = "HTO_H3", cex.use = 0.5, col.use = cols)
GenePlot(hashedCells, gene1 = "HTO_H2", gene2 = "HTO_H3", cex.use = 0.5, col.use = cols)

dev.off()

hashedCells <- SetAllIdent(pbmc_hashtag,"hto_classification")


