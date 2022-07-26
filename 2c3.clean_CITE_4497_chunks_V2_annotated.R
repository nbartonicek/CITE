library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
#library(DropletUtils)
library(scImpute)



dims=40
resolution=0.8
perplexity=20
conditions<-paste(dims,resolution,perplexity,sep="_")


homeDir="/share/ScratchGeneral/nenbar/projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4497-1_CHUNKS2"
sampleName="4497-1_CHUNKS2_CITE"
treatment="v2"
#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))
rawDir=paste0(homeDir,"raw_files/",projectname,"/",sampleName,"/raw_data/")
system(paste0("mkdir -p ",rawDir))




#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))

cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,".tsv"), 
                       sep = ",", header = TRUE, row.names = 1)
cite <- cite[!(row.names(cite) %in% c("bad_struct","no_match","total_reads")),]
rownames(cite) <- paste0("CITE_", rownames(cite))
cite <-cite[,order(colnames(cite))]



tSNEFile=paste0(homeDir,"project_results/Robjects/annotated_4497-1_chunks2.RData")
SeuratX <- readRDS(file = tSNEFile)



#save the raw data as a data frame
raw_data=as.data.frame(as.matrix(citeCells@raw.data))
filtered_data=raw_data %>% select(colnames(citeCells@data))
write.csv(filtered_data,paste0(rawDir,"raw_counts.csv"),quote=F)

#save the data in the mtx format
sparse.gbm <- Matrix(as.matrix(filtered_data) , sparse = T )
writeMM(obj = sparse.gbm, file=paste0(rawDir,"matrix.mtx"))

# save genes and cells names
write(x = rownames(filtered_data), file = paste0(rawDir,"genes.tsv"))
write(x = colnames(filtered_data), file = paste0(rawDir,"barcodes.tsv"))

imputedDir=paste0(homeDir,"raw_files/",projectname,"/",sampleName,"/imputed_data/")
system(paste0("mkdir -p ",imputedDir))

labels<-citeCells@meta.data$res.0.6
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


treatment="v2"


earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_",treatment,"_earlySeurat.Rdata")
seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(cite) ]


  #combine the mouse and human
  combined<-CollapseSpeciesExpressionMatrix(seuratMatClean,prefix.1="GRCh38_",prefix.controls="mm10__",features.controls.toKeep=500)
  
  citeCells <- CreateSeuratObject(
    raw.data = seuratMatClean, 
    min.cells = 3, 
    min.genes = 100, 
    project = sampleName
  )

  save(citeCells,file=earlySeuratFile)

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
  citeCells <- JackStraw(object = citeCells, num.pc = 50,num.replicate = 100)
  citeCells <- JackStrawPlot(object = citeCells, PCs = 1:20)

  save(citeCells,file=jackstrawFile)
}else{load(jackstrawFile)}


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

pdf(paste0(figDir,sampleName,"_",conditions,"_",treatment,"_tsne.pdf"),width=12,height=8)

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

genes<-c("nGene","percent.mito",paste0("GRCh38_",c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,sampleName,"_",conditions,"_feature.pdf"),width=12,height=16)
FeaturePlot(object = citeCells, features.plot = genes, nCol = 4, pt.size = 0.5,cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

dev.off()
#}


#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)
CITEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_CITE.Rdata")

save(citeCells,file=CITEFile)

# visualized on a tSNE generated using the ADT distance matrix.
pdf(paste0(figDir,"all_CITE_",sampleName,".pdf"),width=10,height=50)
FeaturePlot(citeCells, features = row.names(cite), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


#do this on t-cells only
epicells <- SubsetData(citeCells, ident.use = c(6,11,15,2,14,8,4,1,5))

epicells_cite <- RunPCA(epicells, pcs.compute = 45, pc.genes = rownames(cite), assay.type = "CITE", 
    pcs.print = 0)

pdf(paste0(figDir,"epicell_CITE_PCA_",sampleName,".pdf"),width=12,height=4)
PCAPlot(epicells_cite, pt.size = 0.5,cols.use=cols)
dev.off()




adt.data <- GetAssayData(epicells_cite, assay.type = "CITE", slot = "data")
adt.dist <- as.matrix(dist(t(adt.data)))

# Before we recluster the data on ADT levels, we'll stash the RNA cluster
# IDs for later
epicells_cite <- StashIdent(epicells_cite, "rnaClusterID")

# Now, we rerun tSNE using our distance matrix defined only on ADT (protein)
# levels.
epicells_cite <- RunTSNE(epicells_cite, distance.matrix = adt.dist)

# We can also rerun clustering using the same distance matrix. We'll start
# with a very coarse clustering (resolution=0.2)
epicells_cite <- FindClusters(epicells_cite, distance.matrix = adt.dist, print.output = FALSE, 
    resolution = 0.8)

# We can compare the RNA and protein clustering, and use this to annotate
# the protein clustering (we could also of course use FindMarkers)
clustering.table <- table(epicells_cite@ident, epicells_cite@meta.data$rnaClusterID)



tsne_rnaClusters <- TSNEPlot(epicells_cite, do.return = TRUE, group.by = "rnaClusterID", 
    pt.size = 0.5)
tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Clustering based on scRNA-seq") + 
    theme(plot.title = element_text(hjust = 0.5))

tsne_adtClusters <- TSNEPlot(epicells_cite, do.return = TRUE, pt.size = 0.5)
tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on ADT signal") + 
    theme(plot.title = element_text(hjust = 0.5))

# Note: for this comparison, both the RNA and protein clustering are
# visualized on a tSNE generated using the ADT distance matrix.
pdf(paste0(figDir,"epicells_CITE_tsne_comparison_",sampleName,".pdf"),width=10,height=4)
plot_grid(tsne_rnaClusters, tsne_adtClusters, ncol = 2)
dev.off()


pdf(paste0(figDir,"epicells_total_feature_plot_cite_",sampleName,"1.pdf"),width=10,height=50)
FeaturePlot(epicells_cite, features = row.names(cite)[c(1:4,7:22,24:97)], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()




################ use only the annotated cells ############

annotatedShort <- SubsetData(SeuratX, ident.use = c("epithelial_2","epithelial_1","epithelial_5","epithelial_4","epithelial_6"))

epicells <- SubsetData(citeCells, cells.use = gsub("CID44971CHUNKS2_","",colnames(annotatedShort@data)))
annotatedEpiShort <- SubsetData(annotatedShort, cells.use = paste0("CID44971CHUNKS2_",colnames(epicells@data)))

#sepicells_cite@meta.data=annotatedEpiShort@ident$rnaClusterID


cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F")


epicells_cite <- RunPCA(epicells, pcs.compute = 45, pc.genes = rownames(cite), assay.type = "CITE", 
    pcs.print = 0)

pdf(paste0(figDir,"epicell_CITE_PCA_",sampleName,".pdf"),width=12,height=4)
PCAPlot(epicells_cite, pt.size = 0.5,cols.use=cols)
dev.off()




adt.data <- GetAssayData(epicells_cite, assay.type = "CITE", slot = "data")
adt.dist <- as.matrix(dist(t(adt.data)))

# Before we recluster the data on ADT levels, we'll stash the RNA cluster
# IDs for later
epicells_cite <- StashIdent(epicells_cite, "rnaClusterID")
epicells_cite@ident$rnaClusterID<-
# Now, we rerun tSNE using our distance matrix defined only on ADT (protein)
# levels.
epicells_cite <- RunTSNE(epicells_cite, distance.matrix = adt.dist)

# We can also rerun clustering using the same distance matrix. We'll start
# with a very coarse clustering (resolution=0.2)
epicells_cite <- FindClusters(epicells_cite, distance.matrix = adt.dist, print.output = FALSE, 
    resolution = 0.8)

# We can compare the RNA and protein clustering, and use this to annotate
# the protein clustering (we could also of course use FindMarkers)
clustering.table <- table(epicells_cite@ident, epicells_cite@meta.data$rnaClusterID)



tsne_rnaClusters <- TSNEPlot(epicells_cite, do.return = TRUE, group.by = "rnaClusterID", 
    pt.size = 0.5)
tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Clustering based on scRNA-seq") + 
    theme(plot.title = element_text(hjust = 0.5))

tsne_adtClusters <- TSNEPlot(epicells_cite, do.return = TRUE, pt.size = 0.5)
tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on ADT signal") + 
    theme(plot.title = element_text(hjust = 0.5))

dev.off()

# Note: for this comparison, both the RNA and protein clustering are
# visualized on a tSNE generated using the ADT distance matrix.
pdf(paste0(figDir,"epicellsGammy_CITE_tsne_comparison_",sampleName,".pdf"),width=10,height=4)
plot_grid(tsne_rnaClusters, tsne_adtClusters, ncol = 2)
dev.off()


pdf(paste0(figDir,"epicellsGammy_total_feature_plot_cite_",sampleName,"1.pdf"),width=10,height=50)
FeaturePlot(epicells_cite, features = row.names(cite)[c(1:4,7:22,24:97)], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()











