library(dplyr)
library(Matrix)
library(devtools)
library(scImpute)
library(data.table)
library(R.utils)
args <- R.utils::commandArgs(asValues=TRUE)
if (!is.null(args[["order"]])){order = as.numeric(args$order)} 

#inDir<-"/share/ScratchGeneral/sunwu/projects/Stromal_heterogeneity_paper/CCA_alignment_cutoff2/Output/Raw_matrices_and_metadata_for_Nenad_imputation/"
#dirs<-list.files(inDir)
#dir <- dirs[order]
#cat(dir)
#cat("\n")
#csv<-list.files(paste0(inDir,dir),full.names=T,pattern="raw_count_expression_matrix")
#annotation<-read.table(paste0(inDir,dir,"/",dir,"_cluster.txt"),sep="\t",header=T,stringsAsFactors=F)
#annotation<-annotation[order(annotation$cell_barcode),]
#outDir<-paste0("/share/ScratchGeneral/nenbar/projects/CITE/raw_files/Sunny_imputed/",dir,"/")
#
#scimpute(csv, 
#       infile = "csv", 
#       outfile = "csv", 
#       out_dir = outDir,
#       drop_thre = 0.5,
#       ncores = 10,
#       labeled = T,
#       labels = annotation$cell_type)


inDir="/share/ScratchGeneral/sunwu/projects/Stromal_heterogeneity_paper/CCA_alignment_cutoff2/Output/Raw_matrices_and_metadata_for_Nenad_imputation/INTEGRATED/"
csv<-list.files(inDir,full.names=T,pattern="raw_count_expression_matrix")
annotation<-read.table(paste0(inDir,"INTEGRATED_cluster.csv"),sep=",",header=T,stringsAsFactors=F)
annotation<-annotation[order(annotation$cell_barcode),]
outDir<-paste0("/share/ScratchGeneral/nenbar/projects/CITE/raw_files/Sunny_imputed/INTEGRATED/")


scimpute(csv, 
       infile = "csv", 
       outfile = "csv", 
       out_dir = outDir,
       drop_thre = 0.5,
       ncores = 6,
       labeled = T,
       labels = annotation$cell_type)
