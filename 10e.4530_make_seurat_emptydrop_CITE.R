library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4530"
sampleName="4530"
treatments=c("tumor_emptydrop","normal_emptydrop")
treatment="normal_emptydrop"
resolution=0.8
perplexity=20

#for(treatment in treatments){

  #treatment="mt0.5_nGene7000"
  figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
  system(paste0("mkdir -p ",figDir))
  rObjectsDir=paste0(homeDir,"project_results/Robjects/")
  system(paste0("mkdir -p ",rObjectsDir))


  tag=""
  if(treatment=="normal_emptydrop"){tag="-N"}

  tSNEFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_tSNE.Rdata")
  load(tSNEFile)


  if(treatment=="normal_emptydrop"){tag=""}else{tag="-N"}
  cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/CID",sampleName,tag,"_CITE_emptydrop.tsv"), 
                         sep = ",", header = TRUE, row.names = 1)
    
  
  if(treatment=="normal_emptydrop"){tag=""}else{tag="-N"}
  cite <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/CID",sampleName,tag,"_CITE.tsv"), 
                         sep = ",", header = TRUE, row.names = 1)
    
  rownames(cite) <- paste0("CITE_", rownames(cite))
  cite <- cite[1:(dim(cite)[1]-3),]
  cite <-cite[,order(colnames(cite))]
  cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
  citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
  citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
  citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)




  pdf(paste0(figDir,"feature_plot_cite_",sampleName,".pdf"),width=16,height=20)
  FeaturePlot(citeCells, features.plot = row.names(cite), min.cutoff = "q05", max.cutoff = "q95", 
      nCol = 5, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
  dev.off()


  pdf(paste0(testFigDir,"ridgeplot.pdf"),width=12,height=30)

  RidgePlot(citeCells, features.plot = row.names(cite), cols.use = cols,nCol = 4)

  dev.off()












  pdf(paste0(figDir,"feature_plot_hash_",sampleName,".pdf"),width=8,height=4)
  FeaturePlot(hashedCells, features.plot = c("HTO_1", "HTO_2"), min.cutoff = "q05", max.cutoff = "q95", 
      nCol = 2, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
  dev.off()


  pdf(paste0(figDir,"cross_hash_",sampleName,".pdf"),width=8,height=4)
  GenePlot(hashedCells, gene1 = "HTO_1", gene2 = "HTO_2", cex.use = 0.5, col.use = cols)
  dev.off()

  #perform analysis for the subsamples
  subsampleL<-list()

  subsampleL[["1"]]=table(citeCells@meta.data$hto_classification_global)

  for( i in 1:9/100){
    cat(i)
    cat("\n")
    hash <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/subsamples/CID",sampleName,tag,"_HTO_",i,".emptydrop.tsv"), 
                         sep = ",", header = TRUE, row.names = 1)
    rownames(hash) <- gsub("HTO","HTO_", rownames(hash))
    hash <- hash[1:(dim(hash)[1]-3),]
    joint_bcs <- intersect(colnames(citeCells@data),colnames(hash))
    #hashedCells_sparse <- hashedCells@data[,joint_bcs]
    hash <- as.matrix(hash[,joint_bcs])
    hash <- hash[1:2,]
    citeCells <- SetAssayData(citeCells, assay.type = "HTO", slot = "raw.data", new.data = hash)
    citeCells <- NormalizeData(citeCells, assay.type = "HTO", normalization.method = "genesCLR")
    citeCells <- ScaleData(citeCells, assay.type = "HTO", display.progress = FALSE)
    citeCells <- HTODemux(citeCells,assay.type = "HTO",positive_quantile =  0.99,print.output = FALSE)
    subsampleL[[as.character(i)]]=table(citeCells@meta.data$hto_classification_global)
  }

  df=do.call("rbind",subsampleL)
  df=as.data.frame(df)
  df=df[order(row.names(df)),]
  df
  df$Negative/df$Singlet

#}
