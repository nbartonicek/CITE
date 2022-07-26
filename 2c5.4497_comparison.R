library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(reshape2)
#library(DropletUtils)
library(ggrepel)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4497-1_CHUNKS2"
sampleName="4497-1_CHUNKS2_CITE"



treatment="comparison"
#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))
rawDir=paste0(homeDir,"raw_files/",projectname,"/",sampleName,"/raw_data/")
system(paste0("mkdir -p ",rawDir))


types=c("old","new")


res<-list()
for(type in types){
    inFile<-paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_emptydrop_CITE_",type,".tsv")
    res[[paste0(type)]]<-read.csv(inFile,sep=",",header=T,row.names=1)
}

#plot one vs other

citeOLD <- res[[1]]
citeOLD <-citeOLD[!(row.names(citeOLD) %in% c("bad_struct","no_match","total_reads")),]
rownames(citeOLD) <- paste0("CITE_", rownames(citeOLD))
citeOLD <-citeOLD[order(row.names(citeOLD)),order(colnames(citeOLD))]

#load the new emptydrops
citeNEW <- res[[2]]
citeNEW <- citeNEW[!(row.names(citeNEW) %in% c("bad_struct","no_match","total_reads")),]
rownames(citeNEW) <- paste0("CITE_", rownames(citeNEW))
citeNEW <-citeNEW[order(row.names(citeNEW)),order(colnames(citeNEW))]

citeNEW<-citeNEW[row.names(citeNEW) %in% row.names(citeOLD),]
citeOLD<-citeOLD[row.names(citeOLD) %in% row.names(citeNEW),]


#plot2 library count
pdf(paste0(figDir,"libraryCounts.pdf"),width=14,height=8)

df<-data.frame(ab=row.names(citeNEW),old=rowSums(citeOLD),new=rowSums(citeNEW))
#df<-data.frame(ab=row.names(citeNEW),old=log10(rowSums(citeOLD)),new=log10(rowSums(citeNEW)))
df=df[!grepl("\\.",df$ab),]
dataM<-melt(df)
colnames(dataM)<-c("CITE","type","counts")
p<-ggplot(df,aes(old,new))
p<-p+geom_point()
p<-p+scale_x_log10()+scale_y_log10()
#p<-p+geom_text(aes(label = ab),hjust = 0, size=4, alpha=0.5,nudge_x = 10)
p<-p+geom_text_repel(data=df, size=4,
      segment.colour = "gray40",
      direction     = "y",
      aes(x = old, y = new, 
        label = ab))
p
dev.off()


#load the old CITE file
treatment="comparison"


earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_",treatment,"_earlySeurat.Rdata")
if(!file.exists(earlySeuratFile)){
  dataH <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/GRCh38")
  dataM <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/mm10")

  seuratH <- Read10X(dataH)
  seuratM <- Read10X(dataM)
  seuratX<-rbind(seuratH,seuratM)
  #seuratMat<-as.matrix(seuratX[,colSums(seuratX)>500])
  #this produces 4153 cells (cutoff 500)

  whitelist=read.table(file=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_whitelist.txt"))
  whitelist=read.table(file=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_emptydrop_whitelist.txt"))

  seuratMat=seuratX[,colnames(seuratX) %in% colnames(citeNEW)]
  #this produces 4754 cells  (FDR 0.01)
  seuratMatClean<-seuratMat
   #add the cite data
  seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(citeOLD) ]


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























#plot1 stats per library

#plot3 counts per antibody PARALLEL
pdf(paste0(figDir,"antibodyCounts.pdf"),width=12,height=8))

dev.off()

#plot4 antibody correlation coefficients PARALLEL - without both-zeroes
pdf(paste0(figDir,"antibodyCounts.pdf"),width=12,height=8))

dev.off()




#4757 cells
emptydrop_whitelist<-read.table("../raw_files/4497-1_CHUNKS2/CITE_fastqs/4497-1_CHUNKS2_emptydrop_whitelist.txt")

#store the old CITE
oldCite=citeCells@assay$CITE@raw.data
#3791 cells, 98 antibodies, 6.35 million reads

#set the new CITE to it and plot side-by-side comparison
cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_emptydrop_CITE.tsv"), 
                       sep = ",", header = TRUE, row.names = 1)
rownames(cite) <- paste0("CITE_", rownames(cite))
cite <- cite[1:(dim(cite)[1]-3),]
cite <-cite[,order(colnames(cite))]


#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)
CITEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_CITE.Rdata")



#tSNE and correlation


cite <- cite[1:(dim(cite)[1]-3),]















cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_emptydrop_CITE.tsv"), 
                       sep = ",", header = TRUE, row.names = 1)
rownames(cite) <- paste0("CITE_", rownames(cite))
cite <- cite[1:(dim(cite)[1]-3),]
cite <-cite[,order(colnames(cite))]


earlySeuratFile=paste0("../project_results/Robjects/",sampleName,"_",treatment,"_earlySeurat.Rdata")
if(!file.exists(earlySeuratFile)){
  dataH <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/GRCh38")
  dataM <- paste0(homeDir,"raw_files/",projectname,"/outs/raw_gene_bc_matrices/mm10")

  seuratH <- Read10X(dataH)
  seuratM <- Read10X(dataM)
  seuratX<-rbind(seuratH,seuratM)
  #seuratMat<-as.matrix(seuratX[,colSums(seuratX)>500])
  #this produces 4153 cells (cutoff 500)

  whitelist=read.table(file=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_whitelist.txt"))
  whitelist=read.table(file=paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_emptydrop_whitelist.txt"))

  seuratMat=seuratX[,colnames(seuratX) %in% colnames(cite)]
  #this produces 4754 cells  (FDR 0.01)

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
  combined<-CollapseSpeciesExpressionMatrix(seuratMatClean,prefix.1="GRCh38_",prefix.controls="mm10__",features.controls.toKeep=500)
  
  citeCells <- CreateSeuratObject(
    raw.data = combined, 
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
  citeCells <- FindClusters(citeCells, dims.use = pcs, resolution = resolution,print.output = FALSE)
  citeCells <- RunTSNE(citeCells, dims.use = pcs, perplexity=perplexity)
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
cite <-citeOLD[,colnames(citeOLD) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
#citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)
CITEFileOLD=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_CITEOLD.Rdata")
save(citeCells,file=CITEFileOLD)


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"1.OLD.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[1:20], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"2.OLD.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[21:40], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()



#add CITE data
cite <-citeNEW[,colnames(citeNEW) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
#citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)
CITEFileNEW=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_CITENEW.Rdata")
save(citeCells,file=CITEFileNEW)


pdf(paste0(figDir,"feature_plot_cite_",sampleName,"1.NEW.pdf"),width=12,height=10)
FeaturePlot(citeCells, features = row.names(cite)[1:20], min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()



pdf(paste0(figDir,"feature_plot_cite_",sampleName,"2.NEW.pdf"),width=12,height=10)
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

RidgePlot(citeCells, features = row.names(cite), cols = cols,nCol = 3)

dev.off()




###### do the clustering based on CITE #####



#save the raw data as a data frame
as.data.frame(as.matrix(citeCells@raw.data))
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









