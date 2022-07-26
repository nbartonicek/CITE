library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)

dims=40
resolution=0.8
perplexities=c(5,10,15,20,50)
perplexity=20


conditions<-paste(dims,resolution,perplexity,sep="_")
cat(conditions)

homeDir="/share/ScratchGeneral/nenbar/projects/CITE/"
homeDir="../../../projects/CITE/"

figDir=paste0(homeDir,"project_results/figures/")
projectname="4515"
sampleName="4515_CITE"
treatment="umap"
#treatment="mt0.5_nGene7000"
testFigDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")

system(paste0("mkdir -p ",testFigDir))

cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/4515_CITE.tsv"), 
                       sep = ",", header = TRUE, row.names = 1)
rownames(cite) <- paste0("CITE_", rownames(cite))
cite <- cite[1:(dim(cite)[1]-3),]
cite <-cite[,order(colnames(cite))]



earlySeuratFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_earlySeurat.Rdata")
load(earlySeuratFile)

tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_tSNE.Rdata")
load(tSNEFile)

pdf(paste0(testFigDir,sampleName,"_",conditions,"_tsne.pdf"),width=12,height=8)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F","red","blue")
TSNEPlot(object = citeCells, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6,
         colors.use = cols
)
dev.off()

citeCells <- RunUMAP(citeCells, reduction.use = "pca", dims.use = 1:100)
pdf(paste0(testFigDir,sampleName,"_umap.pdf"),width=12,height=8)
DimPlot(citeCells, reduction.use = "umap",colors.use = cols)
dev.off()


genes<-c("nGene","percent.mito","ACTB","EPCAM","KRT8","KRT5","COL1A1","PDGFRA","ACTA2","CD34","PECAM1","PTPRC","CD68","FCGR3A","CD19","CD3D","FOXP3")
pdf(paste0(testFigDir,sampleName,"_feature.pdf"),width=12,height=20)
FeaturePlot(object = citeCells, features.plot = genes, nCol = 3, pt.size = 0.5,cols.use = c("grey", "blue"), 
            reduction.use = "umap")

dev.off()
#}

#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)


pdf(paste0(testFigDir,"feature_plot_cite_",sampleName,".pdf"),width=12,height=20)
FeaturePlot(citeCells, features.plot = row.names(cite), reduction.use = "umap", min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(testFigDir,"ridgeplot.pdf"),width=12,height=30)

RidgePlot(citeCells, features.plot = row.names(cite), cols.use = cols,nCol = 4,reduction.use = "umap")

dev.off()



genes<-c("nGene","percent.mito","ACTB","EPCAM","KRT8","KRT5","COL1A1","PDGFRA","PDGFRA","ACTA2","CD34","PECAM1","PTPRC","CD68","FCGR3A","CD19","CD3D","FOXP3")
pdf(paste0(testFigDir,sampleName,"_",conditions,"_feature.pdf"),width=12,height=20)
FeaturePlot(object = citeCells, features.plot = genes, nCol = 3, pt.size = 0.5,cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

dev.off()
#}

