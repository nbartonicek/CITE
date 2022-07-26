library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
#library(DropletUtils)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4290B"
sampleName="4290B_HTO"
treatment="v2"

dims=40
resolution=0.8
perplexity=20
conditions<-paste(dims,resolution,perplexity,sep="_")



#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))

hash <- read.csv(paste0(homeDir,"raw_files/",projectname,"/HTO_fastqs/",projectname,"_all.tsv"), 
                       sep = ",", header = TRUE, row.names = 1)

hash <- hash[!(row.names(hash) %in% c("bad_struct","no_match","total_reads")),]

#rownames(hash) <- paste0("HTO_", rownames(hash))
hash <-hash[,order(colnames(hash))]
hash <-hash[grep("HTO",row.names(hash)),]



cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_CITE.tsv"), 
                       sep = ",", header = TRUE, row.names = 1)

cite <- cite[!(row.names(cite) %in% c("bad_struct","no_match","total_reads")),]

rownames(cite) <- paste0("CITE_", rownames(cite))
cite <-cite[!grepl("\\.",row.names(cite)),]
cite <-cite[,order(colnames(cite))]




earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_",treatment,"_earlySeurat.Rdata")
if(!file.exists(earlySeuratFile)){
  
  dataH <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/GRCh38")
  #dataHM <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/GRCh38_mm10")
  dataM <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/mm10")

  seuratH <- Read10X(dataH)
  #seuratHM <- Read10X(dataHM)
  seuratM <- Read10X(dataM)
  seuratX<-rbind(seuratH,seuratM)
  #seuratMat<-as.matrix(seuratX[,colSums(seuratX)>500])
  #this produces 4153 cells (cutoff 500)

  whitelist=read.table(file=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_whitelist.txt"))

  seuratMat=seuratX[,colnames(seuratX) %in% whitelist$V1]

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
  #humanCells<-colnames(seuratMat[,speciesRatio>=0.5])
  #
  ##species list
  #speciesList<-speciesRatio[speciesRatio>=0.9|speciesRatio<=0.1]
  #speciesNames<-ifelse(speciesList>0.5,"human","mouse")
  #speciesColor<-ifelse(speciesList>0.5,"darkblue","red")
  #save(speciesColor,file="speciesColor.Rdata")
  #save(speciesRatio,file="speciesRatio.Rdata")

  #add the cite data
  seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(hash) ]
  seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(cite) ]

  #mm10MT<-seuratX[row.names(seuratX)=="mm10___mt-Nd1",]
  #mm10MT=mm10MT[names(mm10MT) %in% colnames(citeCells@data)]
  #citeCells@meta.data$percent.mouse.mito=mm10MT
  #combine the mouse and human
  combined<-CollapseSpeciesExpressionMatrix(seuratMatClean,prefix.1="GRCh38_",prefix.controls="mm10__",features.controls.toKeep=500)

  citeCells <- CreateSeuratObject(
    raw.data = combined, 
    min.cells = 3, 
    min.genes = 100, 
    project = sampleName
  )

  save(citeCells,file=earlySeuratFile)
} else {load(earlySeuratFile)}


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
  citeCells <- JackStraw(object = citeCells, num.pc = 100,num.replicate = 100)
  citeCells <- JackStrawPlot(object = citeCells, PCs = 1:20)

  save(citeCells,file=jackstrawFile)
}else{load(jackstrawFile)}
pdf(paste0(figDir,"vlnplot_",sampleName,".pdf"),width=12,height=6)
VlnPlot(object = citeCells, features = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

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

genes<-c("nGene","percent.mito",paste0(c("ACTB","HLA-A","ATP1B3","CD47","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,sampleName,"_",treatment,"_feature.pdf"),width=12,height=16)
FeaturePlot(object = citeCells, features.plot = genes, nCol = 4, pt.size = 0.5,cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

dev.off()
#}


#pdf(paste0(figDir,sampleName,"_tSNE_with_predefined_genes_human.pdf"),width=8,height=8)
pdf(paste0(figDir,"featurePlot_",sampleName,"_",treatment,".pdf"),width=12,height=8)
FeaturePlot(object = citeCells, pt.size = 0.2,    features = genes, cols = c("grey", "blue"), 
            reduction = "tsne",label.size = 0.5)
dev.off()

citeList<-list()
citeCells <- SetAllIdent(citeCells,"res.0.8")
hash <- hash[1:2,]

clusters<-unique(citeCells@meta.data$res.0.8)
clusters<-clusters[clusters!="9"]



for(cluster in clusters){
  cat(cluster)
  cat("\n")
  tempCite<-SubsetData(object = citeCells, ident.use=cluster)
  hashTemp <-hash[,colnames(hash) %in% colnames(tempCite@data) ]
  tempCite <- SetAssayData(tempCite, assay.type = "HTO", slot = "raw.data", new.data = hashTemp)
  tempCite <- NormalizeData(tempCite, assay.type = "HTO", normalization.method = "genesCLR")
  positive_quantile=0.99
  tempCite <- HTODemux(tempCite,assay.type = "HTO",init_centers=5,positive_quantile =  positive_quantile, print.output = TRUE)
  tempCite <- SetAllIdent(tempCite,"hto_classification")
  citeList[[cluster]]<-tempCite
}

citeCellsCombined <- MergeSeurat(unlist(citeList))
combinedCite <- Reduce(function(x,y){MergeSeurat(x,y)}, citeList)


citeCells@meta.data$ <- SetAllIdent(citeCells,"hto_classification_global")
cols<-c("#b2182b", "gray", "#053061")
pdf(paste0(figDir,"hto_Combined_classification_global_PQ.",positive_quantile,".",sampleName,".pdf"),width=12,height=8)
TSNEPlot(combinedCite,group.by = "hto_classification_global", 
         pt.size = 1, 
         label.size = 6,
         colors.use = cols)
dev.off()












temp=combinedCite
combinedCite <- SetAllIdent(combinedCite,"hto_classification_global")
table(combinedCite@meta.data$hto_classification_global)


cluster="1"
pdf(paste0(figDir,"ridgeplotCluster",cluster,"_PQ.",positive_quantile,".",sampleName,".pdf"),width=12,height=4)
RidgePlot(citeList[[cluster]],features.plot = rownames(GetAssayData(citeList[[cluster]],assay.type = "HTO"))[1:2],nCol = 2)
dev.off()

cluster="7"
pdf(paste0(figDir,"ridgeplotCluster",cluster,"_PQ.",positive_quantile,".",sampleName,".pdf"),width=12,height=4)
RidgePlot(citeList[[cluster]],features.plot = rownames(GetAssayData(citeList[[cluster]],assay.type = "HTO"))[1:2],nCol = 2)
dev.off()

#plot distributions for HTO1, HTO2

#plot distributions of differences




#add CITE data
hash <-hash[,colnames(hash) %in% colnames(citeCells@data) ]
hash <- hash[1:2,]
citeCells <- SetAssayData(citeCells, assay.type = "HTO", slot = "raw.data", new.data = hash)
citeCells <- NormalizeData(citeCells, assay.type = "HTO", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, vars.to.regress="nUMI",assay.type = "HTO", do.par = T, display.progress = T)
positive_quantile=0.99
citeCells <- HTODemux(citeCells,init_centers=7,assay.type = "HTO",positive_quantile =  positive_quantile, print.output = TRUE)

citeCells <- SetAllIdent(citeCells,"hto_classification_global")

pdf(paste0(figDir,"feature_plot_hash4_PQ.",positive_quantile,".",sampleName,".pdf"),width=8,height=4)
FeaturePlot(citeCells, features.plot = c("HTO1", "HTO2"), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 2, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()

pdf(paste0(figDir,"ridgeplot_PQ.",positive_quantile,".",sampleName,".pdf"),width=12,height=4)
RidgePlot(citeCells,features.plot = rownames(GetAssayData(citeCells,assay.type = "HTO"))[1:2],nCol = 2)
dev.off()

citeCells <- SetAllIdent(citeCells,"hto_classification")
cols<-c("#b2182b", "pink", "purple", "#053061", 
                         "orange", "#1c9099", 
                        "gray")

pdf(paste0(figDir,"tsneplot_PQ.",positive_quantile,".",sampleName,".pdf"),width=12,height=8)
TSNEPlot(citeCells,group.by = "hto_classification", 
         pt.size = 1, 
         label.size = 6,
         colors.use = cols)
dev.off()

citeCells <- SetAllIdent(citeCells,"hto_classification_global")
cols<-c("#b2182b", "gray", "#053061")
pdf(paste0(figDir,"hto_classification_global_PQ.",positive_quantile,".",sampleName,".pdf"),width=12,height=8)
TSNEPlot(citeCells,group.by = "hto_classification_global", 
         pt.size = 1, 
         label.size = 6,
         colors.use = cols)
dev.off()



cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
cite <-cite[!grepl("HTO",row.names(cite)),]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)

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

a=sort(rowSums(cite))
df=data.frame(count=a,marker=gsub("CITE_","",names(a)))

df$marker<-factor(df$marker,levels=gsub("CITE_","",names(a)))
pdf(paste0(figDir,"reads_nonlog.pdf"),width=16,height=4)
p<-ggplot(df,aes(marker,count))+geom_point()
p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
dev.off()


#do hashed cells have equal distribution across all the clusters

negs<-split(citeCells@meta.data$hto_classification_global,citeCells@meta.data$res.0.8)
negPercentage<-lapply(negs,function(x){sum(x=="Negative")/length(x)})

data<-do.call(negPercentage,"rbind")


