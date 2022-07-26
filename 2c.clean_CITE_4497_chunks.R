library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(ggplot2)
library(reshape2)
library(ggpubr)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4497-1_CHUNKS2"
sampleName="4497-1_CHUNKS2_CITE"
treatment="original"
#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))

cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,"_emptydrops.tsv"), 
                       sep = ",", header = TRUE, row.names = 1)
cite <- cite[!(row.names(cite) %in% c("bad_struct","no_match","total_reads")),]
#rownames(cite) <- paste0("CITE_", rownames(cite))
cite <-cite[,order(colnames(cite))]
cite <-cite[-grep("\\.",row.names(cite)),]

earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_earlySeurat.Rdata")
if(!file.exists(earlySeuratFile)){
  dataH <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/GRCh38")
  dataM <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/mm10")

  seuratH <- Read10X(dataH)
  seuratM <- Read10X(dataM)
  seuratX<-rbind(seuratH,seuratM)
  #seuratMat<-as.matrix(seuratX[,colSums(seuratX)>500])
  #this produces 4153 cells (cutoff 500)
  whitelist=read.table(file=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_emptydrop_whitelist.txt"))
  cite <-cite[,colnames(cite) %in% whitelist$V1]

  seuratMat=seuratX[,colnames(seuratX) %in% colnames(cite)]
  #seuratMat=seuratX[,colnames(seuratX) %in% colnames(citeCells@data)]
 
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
  cite <-cite[,colnames(cite) %in% names(speciesNames)]

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

cite$ab<-row.names(cite)
dataM<-melt(cite)
species<-data.frame(bc=names(speciesNames),species=speciesNames)
dataM<-merge(dataM,species,by.x="variable",by.y="bc")
dataM$species<-factor(dataM$species,levels=c("mouse","human"))
specMedian<-dataM[dataM$species=="human",]
specMedianValues<-split(specMedian,specMedian$ab)
abMedians<-sapply(specMedianValues,function(x){median(x$value)})

dataM$ab<-factor(dataM$ab,levels=names(sort(abMedians)))
#dataM$value=log10(dataM$value+1)
pdf(paste0("cite_mouse_human_",sampleName,"2.pdf"),width=14,height=18)
p<-ggplot(dataM,aes(species,value))
p<-p+geom_violin(aes(fill=species),draw_quantiles = c(0.5))
p<-p+ facet_wrap(~ab, ncol=8)
p<-p+ scale_y_log10()
p<-p+scale_fill_manual(values=c("darkblue", "darkred"))
p<-p+stat_compare_means(aes(group = species), label = "p.signif", label.y = 10, label.x = 1.5) 
p
dev.off()



mouseBC<-names(speciesNames)[speciesNames=="mouse"]
humanBC<-names(speciesNames)[speciesNames=="human"]
citeMouse<-cite[,colnames(cite) %in% mouseBC]
citeHuman<-cite[,colnames(cite) %in% humanBC]

citeSpecies<-split(cite)


mito.genes <- grep(pattern = "^MT-|^mm10---mt-", x = rownames(x = citeCells@assays$RNA@data), value = TRUE)
percent.mito <- Matrix::colSums(citeCells@assays$RNA@counts[mito.genes, ])/Matrix::colSums(citeCells@assays$RNA@counts)
citeCells <- AddMetaData(object = citeCells, metadata = percent.mito, col.name = "percent.mito")

pdf(paste0(figDir,"vlnplot_",sampleName,".pdf"),width=12,height=6)
VlnPlot(object = citeCells, features = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

citeCells <- subset(x = citeCells, subset = nFeature_RNA > 250 & percent.mito < 0.25)


pdf(paste0(figDir,"statsS3_",sampleName,".pdf"),width=12,height=4)
par(mfrow = c(1, 3))
plot(citeCells@meta.data$nCount_RNA,citeCells@meta.data$nFeature_RNA,pch=15,cex=0.6)
plot(citeCells@meta.data$nCount_RNA,citeCells@meta.data$percent.mito,pch=15,cex=0.6)
plot(citeCells@meta.data$nFeature_RNA,citeCells@meta.data$percent.mito,pch=15,cex=0.6)
dev.off()

citeCells <- NormalizeData(object = citeCells, normalization.method = "LogNormalize", scale.factor = 10000)
citeCells <- FindVariableFeatures(citeCells, do.plot = T, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
#length(x = VariableFeatures(object = citeCells))
#3995 features... seems a lot

jackstrawFile=paste0("../project_results/Robjects/",sampleName,"_jackstrawS3.Rdata")

citeCells <- ScaleData(citeCells, vars.to.regress=c("nCount_RNA","nFeature_RNA"),display.progress = FALSE)
citeCells <- RunPCA(citeCells, dims = 50,pcs.print = 0,npcs = 100)
#citeCells <- ProjectPCA(object = citeCells, pcs.store = 100, do.print = FALSE)
citeCells <- JackStraw(object = citeCells, dims = 50,num.replicate = 100)
citeCells <- ScoreJackStraw(object = citeCells, dims = 1:50)



pdf(paste0(figDir,"PCA_heatmap_plot_",sampleName,".pdf"),width=16,height=30)
PCASigGenes(object = citeCells, pcs.use = 1:30)
dev.off()

pdf(paste0(figDir,"PCA_jackstraw_plot_",sampleName,".pdf"),width=12,height=40)
JackStrawPlot(object = citeCells, dims = 1:50)
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


