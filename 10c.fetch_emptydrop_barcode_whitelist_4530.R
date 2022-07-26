library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(DropletUtils)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="4530"
sampleName="4530"
treatments=c("tumor","normal")

results<-list()
for(treatment in treatments){
	cat(treatment)
  #treatment="mt0.5_nGene7000"
  figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
  system(paste0("mkdir -p ",figDir))
  rObjectsDir=paste0(homeDir,"project_results/Robjects/")
  system(paste0("mkdir -p ",rObjectsDir))


  tag="-N"
  if(treatment=="normal"){tag=""}

  dataH <- paste0(homeDir,"raw_files/4530/count_data/CID4530",tag,"/GRCh38")
  dataM <- paste0(homeDir,"raw_files/4530/count_data/CID4530",tag,"/mm10")

  seuratH <- read10xCounts(dataH,col.names=T)
  seuratM <- read10xCounts(dataM,col.names=T)
  seuratX<-rbind(counts(seuratH),counts(seuratM))

  e.out<-emptyDrops(seuratX)
  is.cell <- e.out$FDR <= 0.01
  whitelist<-colnames(seuratX)[is.cell]
  whitelist=whitelist[!is.na(whitelist)]
  whitelist=gsub("-.*","",whitelist)
  results[[treatment]]<-table(Limited=e.out$Limited, Significant=is.cell)
  if(treatment=="normal"){tag="-N"}else{tag=""}
  write(whitelist,file=paste0(homeDir,"/raw_files/4530/CITE_fastqs/CID4530",tag,"_emptydrop_whitelist.txt"))
}