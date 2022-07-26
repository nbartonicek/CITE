library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(scImpute)
library(data.table)

homeDir="/share/ScratchGeneral/nenbar/projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4497"
sample="CID4497-PDX_Hash"
treatment="original"
treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sample,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))
rObjectsDir=paste0(homeDir,"project_results/Robjects/")
system(paste0("mkdir -p ",rObjectsDir))
outDir=paste0(homeDir,"project_results/Robjects/doubletFinder/")
system(paste0("mkdir -p ",outDir))
rawDir=paste0(homeDir,"raw_files/",projectname,"/",sample,"/raw_data/")
imputedDir=paste0(homeDir,"raw_files/",projectname,"/",sample,"/imputed_data/")
system(paste0("mkdir -p ",rawDir))
system(paste0("mkdir -p ",imputedDir))

#load the annotated data
load(paste0(homeDir,"project_results/Robjects/Human_HASHSEQ_S1_hd2_F_clean.Rdata"))
load(paste0(rObjectsDir,"hashedCells_with_HTO.Rdata"))
load(paste0(rObjectsDir,"imputedHashedCells_with_HTO.Rdata"))
#4961 filtered from 5279 cells!!!

#save the raw data as a data frame
#raw_data=as.data.frame(hashedCells@raw.data)
#filtered_data=raw_data %>% select(colnames(hashedCells@data))
#write.csv(filtered_data,paste0(rawDir,"raw_counts.csv"),quote=F)
#
##save the data in the mtx format
#sparse.gbm <- Matrix(as.matrix(filtered_data) , sparse = T )
#writeMM(obj = sparse.gbm, file=paste0(rawDir,"matrix.mtx"))
#
## save genes and cells names
#write(x = rownames(filtered_data), file = paste0(rawDir,"genes.tsv"))
#write(x = colnames(filtered_data), file = paste0(rawDir,"barcodes.tsv"))

#save the raw data as a data frame
raw_data=as.data.frame(hashedCells@raw.data)
filtered_data=raw_data %>% select(colnames(hashedCells@data))
write.csv(filtered_data,paste0(rawDir,"raw_counts.csv"),quote=F)

#save the data in the mtx format
sparse.gbm <- Matrix(as.matrix(filtered_data) , sparse = T )
writeMM(obj = sparse.gbm, file=paste0(rawDir,"matrix.mtx"))

# save genes and cells names
write(x = rownames(filtered_data), file = paste0(rawDir,"genes.tsv"))
write(x = colnames(filtered_data), file = paste0(rawDir,"barcodes.tsv"))


labels<-hashedCells@meta.data$res.0.8
scimpute(paste0(rawDir,"raw_counts.csv"), 
         infile = "csv", 
         outfile = "csv", 
         out_dir = imputedDir,
         labeled = T,
         labels=labels,
         drop_thre = 0.5,
         ncores = 10)


dataImp<-read.csv(paste0(imputedDir,"scimpute_count.csv"),sep=",",row.names=1)
seuratMatClean<-sapply(dataImp,as.integer)
row.names(seuratMatClean)<-row.names(dataImp)


#add the HASH data
#seuratMatClean=dataImp
hash <- read.csv(paste0(homeDir,"raw_files/",projectname,"/Human_HASHSEQ_S1_hd2_F.csv"), 
                     sep = ",", header = TRUE, row.names = 1)
rownames(hash) <- paste0("HTO_", rownames(hash))
hash <- hash[1:(dim(hash)[1]-3),]
hash <-hash[colnames(hash) %in% colnames(seuratMatClean) ]
hash <-hash[,order(colnames(hash))]

seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(hash) ]
dataImp<-dataImp[,colnames(dataImp) %in% colnames(hash) ]


if(!file.exists(paste0(rObjectsDir,"imputedHashedCells.Rdata"))){
  #combine the mouse and human
  #combined<-CollapseSpeciesExpressionMatrix(seuratMatClean,prefix.1="GRCh38_",prefix.controls="mm10_",features.controls.toKeep=500)

  hashedCells <- CreateSeuratObject(
    raw.data = seuratMatClean, 
    min.cells = 3, 
    min.genes = 100, 
    project = "4497_hash"
  )


  mito.genes <- grep(pattern = "^MT-|^mm10___mt-", x = rownames(x = hashedCells@data), value = TRUE)
  percent.mito <- Matrix::colSums(hashedCells@raw.data[mito.genes, ])/Matrix::colSums(hashedCells@raw.data)

  hashedCells <- AddMetaData(object = hashedCells, metadata = percent.mito, col.name = "percent.mito")
  hashedCells <- FilterCells(object = hashedCells, subset.names = c("nGene", "nUMI","percent.mito"), low.thresholds = c(200, 500,-Inf), high.thresholds = c(Inf, Inf,0.3))

  #hashedCells <- AddMetaData(object = hashedCells, metadata = speciesNames, col.name = "species")
  pdf(paste0(figDir,"vlnplotImputed_",sample,".pdf"),width=12,height=8)
  VlnPlot(object = hashedCells, 
          features.plot = c("nGene", "nUMI", "percent.mito"), 
          nCol = 3)
  dev.off()
  hashedCells <- NormalizeData(object = hashedCells)

  hashedCells <- FindVariableGenes(hashedCells, do.plot = T, y.cutoff = 0.5)

  hashedCells <- ScaleData(hashedCells, vars.to.regress=c("nUMI","nGene"),display.progress = FALSE)

  hashedCells <- RunPCA(hashedCells, pcs.print = 0)
  pdf(paste0(figDir,"PCA_bowl_plotImputed_",sample,".pdf"),width=12,height=8)
  PCElbowPlot(hashedCells)
  dev.off()
  hashedCells <- FindClusters(hashedCells, dims.use = 1:15, print.output = FALSE)
  save(hashedCells,file=paste0(rObjectsDir,"imputedHashedCellsHalfdone.Rdata"))
  hashedCells <- RunTSNE(hashedCells, dims.use = 1:15)
  pdf(paste0(figDir,"TSNEImputed_",sample,".pdf"),width=12,height=8)

  cols<-c("#b2182b", "#1F618D", "#053061", 
                          "#bebada", "#f4a582", "#1c9099", 
                          "#85929E", "#9B59B6", "#74add1", 
                          "#de77ae", "#76D7C4", "#b8e186", 
                          "#4393c3", "#c51b7d", "#b8e186", 
                          "#17A589", "#A93226", "#F5B041", "#F1C40F")
  hashedCells <- SetAllIdent(hashedCells,id = "res.0.8")

  TSNEPlot(object = hashedCells, 
           do.label = T, 
           pt.size = 0.5, 
           label.size = 6,
           colors.use = cols
  )
  dev.off()

  save(hashedCells,file=paste0(rObjectsDir,"imputedHashedCells.Rdata"))
} else {load(paste0(rObjectsDir,"imputedHashedCells.Rdata"))}
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
save(hashedCells,file=paste0(rObjectsDir,"imputedHashedCells_with_HTO.Rdata"))
pdf(paste0(figDir,"ridgeplot2imputed_",sample,".pdf"),width=12,height=4)
RidgePlot(hashedCells,features.plot = rownames(GetAssayData(hashedCells,assay.type = "HTO"))[1:3],nCol = 3)
dev.off()

pdf(paste0(figDir,"feature_plot_hash2Imputed_",sample,".pdf"),width=12,height=4)
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

cleanData<-hashedCells@assay$HTO@raw.data
cleanData<-rbind(cleanData,hashedCells@assay$HTO@data)

#hashedCells <- SetAllIdent(pbmc_hashtag,"hto_classification")

