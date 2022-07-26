library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
#library(DropletUtils)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="EL34C35"
sampleName="EL34C35_CITE"
treatment="v2"
#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))
rawDir=paste0(homeDir,"raw_files/",projectname,"/",sampleName,"/raw_data/")
system(paste0("mkdir -p ",rawDir))

cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",sampleName,".tsv"), 
                       sep = ",", header = TRUE, row.names = 1)
cite <- cite[!(row.names(cite) %in% c("bad_struct","no_match","total_reads")),]
rownames(cite) <- paste0("CITE_", rownames(cite))
cite <-cite[,order(colnames(cite))]


earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_",treatment,"_earlySeurat.Rdata")
if(!file.exists(earlySeuratFile)){
  
  dataH <- paste0("/directflow/GWCCGPipeline/projects/sequencing/181214_A00152_0073_BH3WH5DRXX/run/count/JOSPOW_AshwinCITE/JOSPOW_AshwinCITE/outs/raw_gene_bc_matrices/GRCh38")
  seuratX <- Read10X(dataH)
  seuratMat<-as.matrix(seuratX[,colSums(seuratX)>500])
  #this produces 4153 cells (cutoff 500)

  #whitelist=read.table(file=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_whitelist.txt"))
  #whitelist=read.table(file=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_emptydrop_whitelist.txt"))

  #seuratMat=seuratX[,colnames(seuratX) %in% whitelist$V1]
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

#citeCells <- AddMetaData(object = citeCells, metadata = speciesNames, col.name = "species")
pdf(paste0(figDir,"vlnplot_",sampleName,".pdf"),width=12,height=8)
VlnPlot(object = citeCells, 
        features.plot = c("nGene", "nUMI", "percent.mito"), 
        nCol = 3)
dev.off()

citeCells <- FilterCells(object = citeCells, subset.names = c("nGene", "percent.mito"), low.thresholds = c(100, -Inf), high.thresholds = c(Inf, 0.25))

#from 4433 to 4244 

jackstrawFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_jackstraw.Rdata")
if(!file.exists(jackstrawFile)){

  #mito.genes <- grep(pattern = "^MT-|^mm10___mt-", x = rownames(x = citeCells@data), value = TRUE)
  #percent.mito <- Matrix::colSums(citeCells@raw.data[mito.genes, ])/Matrix::colSums(citeCells@raw.data)

  #citeCells <- AddMetaData(object = citeCells, metadata = percent.mito, col.name = "percent.mito")
  #citeCells <- FilterCells(object = citeCells, subset.names = c("nGene", "percent.mito","nUMI"), low.thresholds = c(100, -Inf, 500), high.thresholds = c(Inf, 0.25,Inf))
  citeCells <- NormalizeData(object = citeCells)
  citeCells <- FindVariableGenes(citeCells, do.plot = T,  mean.function = ExpMean,dispersion.function = LogVMR, x.low.cutoff = 0.125, x.high.cutoff = 3,y.cutoff = 0.5)
  citeCells <- ScaleData(citeCells, vars.to.regress=c("nUMI","nGene"),display.progress = FALSE)
  citeCells <- RunPCA(citeCells, pcs.print = 0,npcs = 100)
  citeCells <- ProjectPCA(object = citeCells, pcs.store = 100, do.print = FALSE)
  citeCells <- JackStraw(object = citeCells, num.pc = 100,num.replicate = 100)
  citeCells <- JackStrawPlot(object = citeCells, PCs = 1:20)

  save(citeCells,file=jackstrawFile)
}else{load(jackstrawFile)}



pdf(paste0(figDir,"PCA_jackstraw_plot_",sampleName,".pdf"),width=12,height=40)
JackStrawPlot(object = citeCells)
dev.off()

pcsData<-as.data.frame(citeCells@dr$pca@jackstraw@overall.p.values)
pcs<-pcsData[pcsData$Score<0.01,"PC"]
  

pdf(paste0(figDir,"PCA_heatmap_plot_",sampleName,".pdf"),width=16,height=30)
PCHeatmap(object = citeCells, pc.use = pcs, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
dev.off()


dims=40

resolution=0.6
perplexity=20
conditions<-paste(dims,resolution,perplexity,sep="_")



tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE.Rdata")
if(!file.exists(tSNEFile)){


  #for(perplexity in perplexities){
  citeCells <- FindClusters(citeCells, dims.use = 1:20, resolution = resolution,print.output = FALSE)
  citeCells <- RunTSNE(citeCells, dims.use = 1:20, perplexity=perplexity)
  save(citeCells,file=tSNEFile)
} else {load(tSNEFile)}


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




#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]

cite[rowSums>]

citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
#citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)
CITEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_CITE.Rdata")
save(citeCells,file=CITEFile)


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"1.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[1:20], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


