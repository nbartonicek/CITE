library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(R.utils)

dims=40
resolution=0.8
perplexities=c(5,10,15,20,50)
perplexity=20


args <- R.utils::commandArgs(asValues=TRUE)
if (!is.null(args[["dims"]])){dims = as.numeric(args$dims)} 
if (!is.null(args[["resolution"]])){resolution = as.numeric(args$resolution)} 
if (!is.null(args[["perplexity"]])){perplexity = as.numeric(args$perplexity)} 


conditions<-paste(dims,resolution,perplexity,sep="_")
cat(conditions)

homeDir="/share/ScratchGeneral/nenbar/projects/CITE/"
homeDir="../../../projects/CITE/"

figDir=paste0(homeDir,"project_results/figures/")
projectname="4515"
sampleName="4515_CITE"
treatment="test"
#treatment="mt0.5_nGene7000"
testFigDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")

system(paste0("mkdir -p ",testFigDir))

cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/4515_CITE.tsv"), 
                       sep = ",", header = TRUE, row.names = 1)
rownames(cite) <- paste0("CITE_", rownames(cite))
cite <- cite[1:(dim(cite)[1]-3),]
cite <-cite[,order(colnames(cite))]



earlySeuratFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_earlySeurat.Rdata")
if(!file.exists(earlySeuratFile)){
  dataH <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/GRCh38")
  dataM <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/mm10")

  seuratH <- Read10X(dataH)
  #seuratM <- Read10X(dataM)
  #seuratX<-rbind(seuratH,seuratM)
  seuratX <- seuratH
  #seuratMat<-as.matrix(seuratX[,colSums(seuratX)>500])
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
  #
  ##species list
  #speciesList<-speciesRatio[speciesRatio>=0.9|speciesRatio<=0.1]
  #speciesNames<-ifelse(speciesList>0.5,"human","mouse")
  #speciesColor<-ifelse(speciesList>0.5,"darkblue","red")
  #save(speciesColor,file="speciesColor.Rdata")
  #save(speciesRatio,file="speciesRatio.Rdata")

  #add the cite data
  #seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(cite) ]


  #combine the mouse and human
  #combined<-CollapseSpeciesExpressionMatrix(seuratMatClean,prefix.1="GRCh38_",prefix.controls="mm10_",features.controls.toKeep=500)
    seuratX <-seuratX[,colnames(seuratX) %in% colnames(cite) ]

  citeCells <- CreateSeuratObject(
    raw.data = seuratX, 
    min.cells = 3, 
    min.genes = 100, 
    project = sampleName
  )

  save(citeCells,file=earlySeuratFile)
} else {load(earlySeuratFile)}

jackstrawFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_jackstraw.Rdata")
if(!file.exists(jackstrawFile)){

  mito.genes <- grep(pattern = "^MT-|^mm10___mt-", x = rownames(x = citeCells@data), value = TRUE)
  percent.mito <- Matrix::colSums(citeCells@raw.data[mito.genes, ])/Matrix::colSums(citeCells@raw.data)

  citeCells <- AddMetaData(object = citeCells, metadata = percent.mito, col.name = "percent.mito")
  citeCells <- FilterCells(object = citeCells, subset.names = c("nGene", "percent.mito","nUMI"), low.thresholds = c(100, -Inf, 500), high.thresholds = c(Inf, 0.25,Inf))
  citeCells <- NormalizeData(object = citeCells)
  citeCells <- FindVariableGenes(citeCells, do.plot = T, y.cutoff = 0.5)
  citeCells <- ScaleData(citeCells, vars.to.regress=c("nUMI","nGene"),display.progress = FALSE)
  citeCells <- RunPCA(citeCells, pcs.print = 0,npcs = 100)
  citeCells <- ProjectPCA(object = citeCells, pcs.store = 100, do.print = FALSE)
  citeCells <- JackStraw(object = citeCells, num.pc = 100,num.replicate = 100)
  citeCells <- JackStrawPlot(object = citeCells, PCs = 1:100)

  save(citeCells,file=jackstrawFile)
}else{load(jackstrawFile)}

tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_tSNE.Rdata")
if(!file.exists(tSNEFile)){


  pcsData<-as.data.frame(citeCells@dr$pca@jackstraw@overall.p.values)
  pcs<-pcsData[pcsData$Score<0.01,"PC"]
  #for(perplexity in perplexities){
  citeCells <- FindClusters(citeCells, dims.use = pcs, resolution = resolution,print.output = FALSE)
  citeCells <- RunTSNE(citeCells, dims.use = pcs, perplexity=perplexity)
  save(citeCells,file=tSNEFile)
  } else {load(tSNEFile)}

pdf(paste0(figDir,sampleName,"_",conditions,"_tsne.pdf"),width=12,height=8)

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

#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)


pdf(paste0(figDir,"feature_plot_cite_",sampleName,".pdf"),width=12,height=20)
FeaturePlot(citeCells, features.plot = row.names(cite), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"ridgeplot.pdf"),width=12,height=30)

RidgePlot(citeCells, features.plot = row.names(cite), cols.use = cols,nCol = 4)

dev.off()



genes<-c("nGene","percent.mito","ACTB","EPCAM","KRT8","KRT5","COL1A1","PDGFRA","PDGFRA","ACTA2","CD34","PECAM1","PTPRC","CD68","FCGR3A","CD19","CD3D","FOXP3")
pdf(paste0(testFigDir,sampleName,"_",conditions,"_feature.pdf"),width=12,height=20)
FeaturePlot(object = citeCells, features.plot = genes, nCol = 3, pt.size = 0.5,cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

dev.off()
#}

