library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4530"
sampleName="4530"
treatments=c("tumor_emptydrop","normal_emptydrop")
treatment="tumor_emptydrop"
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

  dataH <- paste0("/share/ScratchGeneral/sunwu/Chromium_10X/count_data_2/CID4530",tag,"/count_CID4520",tag,"_GRCh38_mm10/outs/raw_gene_bc_matrices/GRCh38")
tag,"/count_CID4520",tag,"_GRCh38_mm10/outs/raw_gene_bc_matrices/mm10")

  seuratH <- Read10X(dataH)
  seuratM <- Read10X(dataM)
  seuratX<-rbind(seuratH,seuratM)

  if(treatment=="normal_emptydrop"){tag=""}else{tag="-N"}
  whitelist<-read.table(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/CID",sampleName,tag,"_emptydrop_whitelist.txt"))
  seuratX=seuratX[,colnames(seuratX) %in% whitelist$V1]
  seuratMat<-as.matrix(seuratX)

  #species specific
  humanIdx<-grepl("GRCh38_",row.names(seuratMat))
  speciesRatio<-apply(seuratMat,2,function(x){sum(x[humanIdx])/sum(x)})

  #higher proportion of human than mouse
  pdf(paste0(figDir,"species_ratio_",sampleName,".pdf"),width=12,height=8)
  plot(density(speciesRatio))
  dev.off()
  pdf(paste0(figDir,"species_ratio_hist_",sampleName,".pdf"),width=12,height=8)
  plot(hist(speciesRatio))
  dev.off()

  seuratMatCleanObject<-paste0("../project_results/Robjects/",sampleName,"_",treatment,"_seuratMatClean.Rdata")
  hashObject<-paste0("../project_results/Robjects/",sampleName,"_",treatment,"_hash.Rdata")

  if(!file.exists(seuratMatCleanObject)){
    #eliminate doublets (only 6 with cutoff of 0.9)
    seuratMatClean<-seuratMat[,speciesRatio>=0.9|speciesRatio<=0.1]

    #species list
    speciesList<-speciesRatio[speciesRatio>=0.9|speciesRatio<=0.1]
    speciesNames<-ifelse(speciesList>0.5,"human","mouse")
    speciesColor<-ifelse(speciesList>0.5,"darkblue","red")

    #add the HASH data
    #reverse tag
    if(treatment=="normal_emptydrop"){tag=""}else{tag="-N"}
    hash <- read.csv(paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/CID",sampleName,tag,"_HTO_emptydrop.tsv"), 
                         sep = ",", header = TRUE, row.names = 1)
    rownames(hash) <- gsub("HTO","HTO_", rownames(hash))
    hash <- hash[1:(dim(hash)[1]-3),]
    hash <-hash[colnames(hash) %in% colnames(seuratMatClean) ]
    #6931 to 6911 cells
    hash <-hash[,order(colnames(hash))]
    seuratMatClean <-seuratMatClean[,colnames(seuratMatClean) %in% colnames(hash) ]

    save(hash,file=hashObject)
    save(seuratMatClean,file=seuratMatCleanObject)
  }else{
    load(hashObject)
    load(seuratMatCleanObject)
  }

  earlySeuratFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_earlySeurat.Rdata")

  if(!file.exists(earlySeuratFile)){
    #combine the mouse and human
    combined<-CollapseSpeciesExpressionMatrix(seuratMatClean,prefix.1="GRCh38_",prefix.controls="mm10_",features.controls.toKeep=500)
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
    citeCells <- FilterCells(object = citeCells, subset.names = c("nGene", "percent.mito","nUMI"), low.thresholds = c(100, -Inf, 500), high.thresholds = c(Inf, 0.2,Inf))
    citeCells <- NormalizeData(object = citeCells)
    citeCells <- FindVariableGenes(citeCells, do.plot = T, y.cutoff = 0.5)
    citeCells <- ScaleData(citeCells, vars.to.regress=c("nUMI","nGene"),display.progress = FALSE)
    citeCells <- RunPCA(citeCells, pcs.print = 0,pcs.compute = 100)
    citeCells <- ProjectPCA(object = citeCells, pcs.store = 100, do.print = FALSE)
    citeCells <- JackStraw(object = citeCells, num.pc = 100,num.replicate = 100)
    citeCells <- JackStrawPlot(object = citeCells, PCs = 1:100)

    save(citeCells,file=jackstrawFile)
  }else{load(jackstrawFile)}


  pdf(paste0(figDir,"vlnplot_",sampleName,"_",treatment,".pdf"),width=12,height=8)
  VlnPlot(object = citeCells, 
          features.plot = c("nGene", "nUMI", "percent.mito"), 
          nCol = 3)
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

  pdf(paste0(figDir,sampleName,"_",treatment,"_tsne.pdf"),width=12,height=8)

  cols<-c("#b2182b", "#1F618D", "#053061", 
                          "#bebada", "#f4a582", "#1c9099", 
                          "#85929E", "#9B59B6", "#74add1", 
                          "#de77ae", "#76D7C4", "#b8e186", 
                          "#4393c3", "#c51b7d", "#b8e186", 
                          "#17A589", "#A93226", "#F5B041", "#F1C40F","darkblue","darkred")
  TSNEPlot(object = citeCells, 
           do.label = T, 
           pt.size = 2, 
           label.size = 6,
           colors.use = cols
  )

  dev.off()


  genes<-c("nGene","percent.mito","ACTB","EPCAM","KRT8","KRT5","COL1A1","PDGFRA","ACTA2","CD34","PECAM1","PTPRC","CD68","FCGR3A","CD19","CD3D","FOXP3","mm10___mt-Co3")
  pdf(paste0(figDir,sampleName,"_feature.pdf"),width=16,height=8,useDingbats=FALSE)
  FeaturePlot(object = citeCells, features.plot = genes, nCol = 6, pt.size = 0.5,cols.use = c("grey", "blue"), 
              reduction.use = "tsne")

  dev.off()



  #do the hashing
  joint_bcs <- intersect(colnames(citeCells@data),colnames(hash))
  #hashedCells_sparse <- hashedCells@data[,joint_bcs]
  hash <- as.matrix(hash[,joint_bcs])
  hash <- hash[1:2,]
  citeCells <- SetAssayData(citeCells, assay.type = "HTO", slot = "raw.data", new.data = hash)
  citeCells <- NormalizeData(citeCells, assay.type = "HTO", normalization.method = "genesCLR")
  citeCells <- ScaleData(citeCells, assay.type = "HTO", display.progress = FALSE)
  citeCells <- HTODemux(citeCells,assay.type = "HTO",positive_quantile =  0.99,print.output = FALSE)
  print (table(citeCells@meta.data$hto_classification_global))

  hashedCells <- SetAllIdent(citeCells,id = "hash_maxID")
  save(hashedCells,file=paste0(rObjectsDir,sampleName,"_",treatment,"_hashedCells_with_HTO.Rdata"))
  pdf(paste0(figDir,"ridgeplot_",sampleName,".pdf"),width=8,height=4)
  RidgePlot(hashedCells,features.plot = rownames(GetAssayData(hashedCells,assay.type = "HTO"))[1:2],nCol = 2)
  dev.off()

  pdf(paste0(figDir,"feature_plot_hash_",sampleName,".pdf"),width=8,height=4)
  FeaturePlot(hashedCells, features.plot = c("HTO_1", "HTO_2"), min.cutoff = "q05", max.cutoff = "q95", 
      nCol = 2, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
  dev.off()


  pdf(paste0(figDir,"cross_hash_",sampleName,".pdf"),width=8,height=4)
  GenePlot(hashedCells, gene1 = "HTO_1", gene2 = "HTO_2", cex.use = 0.5, col.use = cols)
  dev.off()

  pdf(paste0(figDir,"hash_annotated_",sampleName,".pdf"),width=12,height=8)

  hashedCells <- SetAllIdent(hashedCells,"hto_classification")

 TSNEPlot(object = citeCells, 
           do.label = F, 
           pt.size = 1, 
           label.size = 6,
           colors.use = cols
  )
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
