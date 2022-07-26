library(Seurat)
library(dplyr)
library(Matrix)
#library(devtools)
#library(DropletUtils)
library(plyr)
library(reshape2)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
RobjectsDir=paste0(homeDir,"project_results/Robjects/")

#projectnames=c( "3838","3946","4040" )
projectnames=c( "3946","4040" )
for(projectname in projectnames){
  cat(projectname)
  cat("\n")
  sampleName=paste0(projectname,"_CITE")
  treatment="v3"

  #treatment="mt0.5_nGene7000"
  figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
  system(paste0("mkdir -p ",figDir))




  #/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_*-1*092019/run01_sunwupipeline/output/count/CID$sampleName/count_CID"$sampleName"_GRCh38_mm10/outs/raw_feature_bc_matrix/
  if(projectname=="3838"){
    inDir=paste0("/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_CITE-16092019/run01_sunwupipeline/output/count/CID",projectname,"/count_CID",projectname,"_GRCh38_mm10/outs/filtered_feature_bc_matrix/")
  } else if (projectname=="3946"){
    inDir=paste0("/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_CITE-19092019/run01_sunwupipeline/output/count/CID",projectname,"/count_CID",projectname,"_GRCh38_mm10/outs/filtered_feature_bc_matrix/")
  } else if (projectname=="4040"){
    inDir=paste0("/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_HTO-17092019/run01_sunwupipeline/output/count/CID",projectname,"/count_CID",projectname,"_GRCh38_mm10/outs/filtered_feature_bc_matrix/")
  }
  if(projectname=="3838"){
    inFile=paste0("/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_CITE-16092019/run01_sunwupipeline/output/seurat/seurat_CID",projectname,"/Output/Rdata/03_seurat_object_processed.RData")
    inFile=paste0(RobjectsDir,"/3838_CITE_v3_tSNE_vst.Rdata")
  } else if (projectname=="3946"){
    inFile=paste0("/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_CITE-19092019/run01_sunwupipeline/output/seurat/seurat_CID",projectname,"/Output/Rdata/03_seurat_object_processed.RData")
    inFile=paste0(RobjectsDir,"/3946_CITE_v3_tSNE_vst.Rdata")
  } else if (projectname=="4040"){
    inFile=paste0("/share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_HTO-17092019/run01_sunwupipeline/output/seurat/seurat_CID",projectname,"/Output/Rdata/03_seurat_object_processed.RData")
    inFile=paste0(RobjectsDir,"/4040_CITE_v3_tSNE_vst.Rdata")
  }
  #cells<-readRDS(inFile)
  load(inFile)
  cells<-citeCells

  matrix_dir = paste0(homeDir,"raw_files/",projectname,"/CITE_fastqs/",projectname,"_CITE.v3/umi_count/")
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  #mat<-mat[colnames(mat) %in% gsub(".*_","",colnames(cells$RNA@data))]
  cite<-as.data.frame(mat)
  cite <-cite[,colnames(cite) %in% gsub(".*_","",colnames(cells$RNA@data)) ]
  rownames(cite) <- paste0("CITE_", rownames(cite))
  rownames(cite) <- gsub("-\\w{15}$","",rownames(cite))
  cite <-cite[,order(colnames(cite))]
  cite <- cite[!(row.names(cite) %in% c("bad_struct","no_match","total_reads")),]
  #colnames(cite)<-paste0("CID",projectname,"_",colnames(cite))
  cite<-cite[,colnames(cite) %in% colnames(cells@assays$RNA@data)]



  cite$id<-row.names(cite)
  allBCs<-colnames(cells@assays$RNA@data)

  citeAll<-as.data.frame(matrix(0,dim(cite)[1],length(allBCs)))
  row.names(citeAll)<-row.names(cite)
  colnames(citeAll)<-allBCs

  citeAll$id<-row.names(citeAll)
  citeL<-list(cite=cite,citeAll=citeAll)

  matrix.df <- ldply(citeL, melt)
  sum.matrix <- acast(matrix.df, id ~ variable, sum)

  citeAll<-sum.matrix
  cite<-cite[,colnames(cite) %in% colnames(cells@assays$RNA@data)]




  citeNames<-gsub("_","-",row.names(citeAll))

  cells[["ADT"]] <- CreateAssayObject(counts = citeAll)
  cells <- NormalizeData(cells, assay = "ADT", normalization.method = "CLR")
  cells <- ScaleData(object = cells, assay = "ADT")

  save(cells,file=paste0(RobjectsDir,projectname,"_v3_CITE.Rdata"))


  #pdf(paste0(figDir,"feature_plot_cite_",sampleName,"1.pdf"),width=12,height=10)
  #FeaturePlot(cells, features = citeNames[1:20], min.cutoff = "q05", max.cutoff = "q95", 
  #    ncol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5,reduction="UMAPC")
  #dev.off()
#
#
  #pdf(paste0(figDir,"feature_plot_cite_",sampleName,"2.pdf"),width=12,height=10)
  #FeaturePlot(cells, features = citeNames[21:40], min.cutoff = "q05", max.cutoff = "q95", 
  #    ncol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5,reduction="UMAPC")
  #dev.off()
#
#
  #pdf(paste0(figDir,"feature_plot_cite_",sampleName,"3.pdf"),width=12,height=10)
  #FeaturePlot(cells, features = citeNames[41:60], min.cutoff = "q05", max.cutoff = "q95", 
  #    ncol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5,reduction="UMAPC")
  #dev.off()
#
#
  #pdf(paste0(figDir,"feature_plot_cite_",sampleName,"4.pdf"),width=12,height=10)
  #FeaturePlot(cells, features = citeNames[61:80], min.cutoff = "q05", max.cutoff = "q95", 
  #    ncol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5,reduction="UMAPC")
  #dev.off()
#
#
  #pdf(paste0(figDir,"feature_plot_cite_",sampleName,"5.pdf"),width=12,height=10)
  #FeaturePlot(cells, features = citeNames[81:100], min.cutoff = "q05", max.cutoff = "q95", 
  #    ncol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5,reduction="UMAPC")
  #dev.off()
#
  #pdf(paste0(figDir,"feature_plot_cite_",sampleName,"5.pdf"),width=12,height=10)
  #FeaturePlot(cells, features = citeNames[101:118], min.cutoff = "q05", max.cutoff = "q95", 
  #    ncol = 4, cols = c("lightgrey", "blue"), pt.size = 0.5,reduction="UMAPC")
  #dev.off()
#
  #adt.markers <- FindAllMarkers(cells, assay = "ADT", only.pos = TRUE)
  #pdf(paste0(figDir,"CITE_garnett_call_ext_major_clusters.pdf"),width=12,height=20)
  #DoHeatmap(cells, features = unique(adt.markers$gene),raster=F, assay = "ADT", angle = 90,group.by="garnett_seurat_cluster_call_subset_PC_C_res.1.2") + NoLegend()
  #dev.off()


}

