library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
#library(DropletUtils)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4497-2_CELLSUS"
sampleName="4497-2_CELLSUS_CITE"
treatment="v2"

dims=40

resolution=1
perplexity=20
conditions<-paste(dims,resolution,perplexity,sep="_")



#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))

cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,"_emptydrops.tsv"), 
                       sep = ",", header = TRUE, row.names = 1)

cite <- cite[!(row.names(cite) %in% c("bad_struct","no_match","total_reads")),]

rownames(cite) <- paste0("GRCh38_CITE_", rownames(cite))
cite <-cite[,order(colnames(cite))]


earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_",treatment,"_earlySeurat_testJoiningCite.Rdata")
if(!file.exists(earlySeuratFile)){
  dataH <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/GRCh38")
  dataM <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/mm10")

  seuratH <- Read10X(dataH)
  seuratM <- Read10X(dataM)
  seuratX<-rbind(seuratH,seuratM)
  #seuratMat<-as.matrix(seuratX[,colSums(seuratX)>500])
  #this produces 4153 cells (cutoff 500)

  whitelist=read.table(file=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_emptydrop_whitelist.txt"))

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
  cite <-cite[,colnames(cite) %in% colnames(seuratMatClean) ]

  temp=rbind(as.matrix(seuratMatClean),as.matrix(cite))
  #combine the mouse and human
  combined<-CollapseSpeciesExpressionMatrix(temp,prefix.1="GRCh38_",prefix.controls="mm10__",features.controls.toKeep=500)
  
  citeCells <- CreateSeuratObject(
    raw.data = combined, 
    min.cells = 3, 
    min.genes = 100, 
    project = sampleName
  )

  save(citeCells,file=earlySeuratFile)
} else {load(earlySeuratFile)}


jackstrawFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_jackstraw_testJoiningCite.Rdata")
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


tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE_testJoiningCite.Rdata")
if(!file.exists(tSNEFile)){


  pcsData<-as.data.frame(citeCells@dr$pca@jackstraw@overall.p.values)
  pcs<-pcsData[pcsData$Score<0.01,"PC"]
  #for(perplexity in perplexities){
  citeCells <- FindClusters(citeCells, dims.use = pcs, resolution = resolution,print.output = FALSE)
  citeCells <- RunTSNE(citeCells, dims.use = pcs, perplexity=perplexity)
  save(citeCells,file=tSNEFile)
} else {load(tSNEFile)}

dev.off()

pdf(paste0(figDir,sampleName,"_",conditions,"_",treatment,"_tsne_testJoiningCite.pdf"),width=12,height=8)

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
pdf(paste0(figDir,sampleName,"_",conditions,"_PCAheatmaps_testJoiningCite.pdf"),width=14,height=30)
PCHeatmap(object = citeCells, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, 
    label.columns = FALSE, use.full = FALSE)
dev.off()

genes<-c("nGene","percent.mito",paste0(c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2","FOXP3","IL2RA")))
pdf(paste0(figDir,sampleName,"_",conditions,"_feature_testJoiningCite.pdf"),width=12,height=16)
FeaturePlot(object = citeCells, features.plot = genes, nCol = 4, pt.size = 0.5,cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

dev.off()
#}


#pdf(paste0(figDir,sampleName,"_tSNE_with_predefined_genes_human.pdf"),width=8,height=8)
#pdf(paste0(figDir,"featurePlot_",sampleName,"_",treatment,"_feature_testJoiningCite.pdf"),width=12,height=8)
#FeaturePlot(object = citeCells, pt.size = 0.2,    features = genes, cols = c("grey", "blue"), 
#            reduction = "tsne")
#dev.off()


#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)
CITEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_CITE_testJoiningCite.Rdata")

save(citeCells,file=CITEFile)


# visualized on a tSNE generated using the ADT distance matrix.
pdf(paste0(figDir,"all_CITE_",sampleName,"_testJoiningCite.pdf"),width=10,height=50)
FeaturePlot(citeCells, features = row.names(cite), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()



pdf(paste0(figDir,"feature_plot_cite_",sampleName,"1.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[1:20], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"2.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[21:40], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"3.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[41:60], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"4.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[61:80], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"5.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[81:93], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"ridgeplot_",sampleName,".pdf"),width=12,height=4)

RidgePlot(citeCells, features = row.names(cite), cols = cols,ncol = 3)

dev.off()




###### do the clustering based on CITE #####

tcells <- SubsetData(citeCells, ident.use = c(0,2,3,4))
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F")

citeCells_cite <- RunPCA(citeCells, pcs.compute = 45, pc.genes = rownames(cite), assay.type = "CITE", 
    pcs.print = 0)

pdf(paste0(figDir,"CITE_PCA_",sampleName,".pdf"),width=12,height=4)
PCAPlot(citeCells_cite, pt.size = 0.5,cols.use=cols)
dev.off()




adt.data <- GetAssayData(citeCells_cite, assay.type = "CITE", slot = "data")
adt.dist <- as.matrix(dist(t(adt.data)))

# Before we recluster the data on ADT levels, we'll stash the RNA cluster
# IDs for later
citeCells_cite <- StashIdent(citeCells_cite, "rnaClusterID")

# Now, we rerun tSNE using our distance matrix defined only on ADT (protein)
# levels.
citeCells_cite <- RunTSNE(citeCells_cite, distance.matrix = adt.dist)

# We can also rerun clustering using the same distance matrix. We'll start
# with a very coarse clustering (resolution=0.2)
citeCells_cite <- FindClusters(citeCells_cite, distance.matrix = adt.dist, print.output = FALSE, 
    resolution = 0.6)

# We can compare the RNA and protein clustering, and use this to annotate
# the protein clustering (we could also of course use FindMarkers)
clustering.table <- table(citeCells_cite@ident, citeCells_cite@meta.data$rnaClusterID)


#current.cluster.ids <- 0:10
## Note, for simplicity we are merging two CD14+ Mono clusters (that differ
## in the expression of HLA-DR genes), and two NK clusters (that differ in
## cell cycle stage)
#new.cluster.ids <- c("CD4 T", "CD14+ Mono", "NK", "B", "CD8 T", "CD34+", "Unknown1", 
#    "CD16+ Mono", "Unknown2", "pDC", "Unknown3")
#cbmc_cite@ident <- plyr::mapvalues(x = cbmc_cite@ident, from = current.cluster.ids, 
#    to = new.cluster.ids)

tsne_rnaClusters <- TSNEPlot(citeCells_cite, do.return = TRUE, group.by = "rnaClusterID", 
    pt.size = 0.5)
tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Clustering based on scRNA-seq") + 
    theme(plot.title = element_text(hjust = 0.5))

tsne_adtClusters <- TSNEPlot(citeCells_cite, do.return = TRUE, pt.size = 0.5)
tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on ADT signal") + 
    theme(plot.title = element_text(hjust = 0.5))

# Note: for this comparison, both the RNA and protein clustering are
# visualized on a tSNE generated using the ADT distance matrix.
pdf(paste0(figDir,"CITE_tsne_comparison_",sampleName,".pdf"),width=12,height=4)
plot_grid(tsne_rnaClusters, tsne_adtClusters, ncol = 2)
dev.off()

pdf(paste0(figDir,"CITE_tsne_",sampleName,".pdf"),width=12,height=4)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F")

TSNEPlot(object = citeCells_cite, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         colors.use = cols
)
dev.off()



#do this on t-cells only

tcells_cite <- RunPCA(tcells, pcs.compute = 45, pc.genes = rownames(cite), assay.type = "CITE", 
    pcs.print = 0)

pdf(paste0(figDir,"tcell_CITE_PCA_",sampleName,".pdf"),width=12,height=4)
PCAPlot(tcells_cite, pt.size = 0.5,cols.use=cols)
dev.off()




adt.data <- GetAssayData(tcells_cite, assay.type = "CITE", slot = "data")
adt.dist <- as.matrix(dist(t(adt.data)))

# Before we recluster the data on ADT levels, we'll stash the RNA cluster
# IDs for later
tcells_cite <- StashIdent(tcells_cite, "rnaClusterID")

# Now, we rerun tSNE using our distance matrix defined only on ADT (protein)
# levels.
tcells_cite <- RunTSNE(tcells_cite, distance.matrix = adt.dist)

# We can also rerun clustering using the same distance matrix. We'll start
# with a very coarse clustering (resolution=0.2)
tcells_cite <- FindClusters(tcells_cite, distance.matrix = adt.dist, print.output = FALSE, 
    resolution = 0.8)

# We can compare the RNA and protein clustering, and use this to annotate
# the protein clustering (we could also of course use FindMarkers)
clustering.table <- table(tcells_cite@ident, tcells_cite@meta.data$rnaClusterID)



tsne_rnaClusters <- TSNEPlot(tcells_cite, do.return = TRUE, group.by = "rnaClusterID", 
    pt.size = 0.5)
tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Clustering based on scRNA-seq") + 
    theme(plot.title = element_text(hjust = 0.5))

tsne_adtClusters <- TSNEPlot(tcells_cite, do.return = TRUE, pt.size = 0.5)
tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on ADT signal") + 
    theme(plot.title = element_text(hjust = 0.5))



pdf(paste0(figDir,"CITE_tsne_",sampleName,".pdf"),width=12,height=4)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F")

TSNEPlot(object = citeCells_cite, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         colors.use = cols
)
dev.off()





# Note: for this comparison, both the RNA and protein clustering are
# visualized on a tSNE generated using the ADT distance matrix.
pdf(paste0(figDir,"tcells_CITE_tsne_comparison_",sampleName,".pdf"),width=10,height=4)
plot_grid(tsne_rnaClusters, tsne_adtClusters, ncol = 2)
dev.off()


pdf(paste0(figDir,"tcells_total_feature_plot_cite_",sampleName,"1.pdf"),width=10,height=50)
FeaturePlot(tcells_cite, features = row.names(cite), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()

pdf(paste0(figDir,"tcells_feature_plot_cite_",sampleName,"1.pdf"),width=10,height=10)
FeaturePlot(tcells_cite, features = row.names(cite)[1:20], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"tcells_feature_plot_cite_",sampleName,"2.pdf"),width=10,height=10)
FeaturePlot(tcells_cite, features = row.names(cite)[21:40], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"tcells_feature_plot_cite_",sampleName,"3.pdf"),width=10,height=10)
FeaturePlot(tcells_cite, features = row.names(cite)[41:60], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"tcells_feature_plot_cite_",sampleName,"4.pdf"),width=10,height=10)
FeaturePlot(tcells_cite, features = row.names(cite)[61:80], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"tcells_feature_plot_cite_",sampleName,"5.pdf"),width=10,height=10)
FeaturePlot(tcells_cite, features = row.names(cite)[81:93], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()



#eliminate the empty cluster
tcells_cite_short <- SubsetData(tcells_cite, ident.use = c(0,1,2,3,4,5,6))


adt.data <- GetAssayData(tcells_cite_short, assay.type = "CITE", slot = "data")
adt.dist <- as.matrix(dist(t(adt.data)))

# Before we recluster the data on ADT levels, we'll stash the RNA cluster
# IDs for later
tcells_cite_short <- StashIdent(tcells_cite_short, "rnaClusterID")

# Now, we rerun tSNE using our distance matrix defined only on ADT (protein)
# levels.
tcells_cite_short <- RunTSNE(tcells_cite_short, distance.matrix = adt.dist)

# We can also rerun clustering using the same distance matrix. We'll start
# with a very coarse clustering (resolution=0.2)
tcells_cite_short <- FindClusters(tcells_cite_short, distance.matrix = adt.dist, print.output = FALSE, 
    resolution = 1)

# We can compare the RNA and protein clustering, and use this to annotate
# the protein clustering (we could also of course use FindMarkers)
clustering.table <- table(tcells_cite_short@ident, tcells_cite_short@meta.data$rnaClusterID)



tsne_rnaClusters <- TSNEPlot(tcells_cite_short, do.return = TRUE, group.by = "rnaClusterID", 
    pt.size = 0.5)
tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Clustering based on scRNA-seq") + 
    theme(plot.title = element_text(hjust = 0.5))

tsne_adtClusters <- TSNEPlot(tcells_cite_short, do.return = TRUE, pt.size = 0.5)
tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on ADT signal") + 
    theme(plot.title = element_text(hjust = 0.5))



# Note: for this comparison, both the RNA and protein clustering are
# visualized on a tSNE generated using the ADT distance matrix.
pdf(paste0(figDir,"tcellsShort_CITE_tsne_comparison_",sampleName,".pdf"),width=10,height=4)
plot_grid(tsne_rnaClusters, tsne_adtClusters, ncol = 2)
dev.off()


pdf(paste0(figDir,"tcellsShort_total_feature_plot_cite_",sampleName,"1.pdf"),width=10,height=50)
FeaturePlot(tcells_cite_short, features = row.names(cite), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()






tcells <- SubsetData(citeCells, ident.use = c(0,2,8))


pdf(paste0(figDir,"t-cellsOnly_",sampleName,"_",conditions,"_",treatment,"_tsne.pdf"),width=8,height=6)

cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F")
TSNEPlot(object = tcells, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         colors.use = cols
)

dev.off()






















