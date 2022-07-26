library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
#library(DropletUtils)

homeDir="../../../projects/CITE/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="MIA"
sampleName="MIA"
treatment="v3"

dims=40
resolution=0.8
perplexity=20
conditions<-paste(dims,resolution,perplexity,sep="_")



#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))

citeCells<-readRDS("/share/ScratchGeneral/ghaale/projects/10x_data/melanoma_project/03_Analysis/RUN_01/Out/Objects/07_Merged__Harmony.rds")
cite<-citeCells@assays$ADT@counts
row.names(cite)[133]<-"C5L2"
cite<-cite[!grepl("unmapped",row.names(cite)),]
citeCells[["ADT"]] <- CreateAssayObject(counts = cite)
citeCells <- NormalizeData(citeCells, assay = "ADT", normalization.method = "CLR")
citeCells <- ScaleData(object = citeCells, assay = "ADT")

  

Idents(object = citeCells) <- "seurat_clusters"

adt.markers <- FindAllMarkers(citeCells, assay = "ADT", only.pos = TRUE)
deABs<-unique(adt.markers$gene)
deABs<-sort(deABs)
pdf(paste0(figDir,"CITE_clusters_merged_harmony.pdf"),width=12,height=20)
p<-DoHeatmap(citeCells, features = unique(adt.markers$gene),raster=F, assay = "ADT", angle = 90) + NoLegend()
print(p)
dev.off()
  

pdf(paste0(figDir,"CITE_clusters_merged_harmony_tsne.pdf"),width=12,height=8)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F","lightgreen","darkblue","lightblue","darkblue","orange","purple","gray","green","black","pink")
DimPlot(object = citeCells, reduction = "umap",
         label = TRUE, 
         pt.size = 2, 
         label.size = 6,
         cols = cols
)
dev.off()



genes<-c("nFeature_RNA","percent.mt",paste0(c("ACTB","CD14","TYR","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,"CITE_clusters_merged_harmony_feature_vst.pdf"),width=12,height=16)
FeaturePlot(object = citeCells, features = genes, ncol = 3, pt.size = 0.5,cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()


#plot all cite data
nABs<-length(deABs)
nPlots<- nABs %/% 20
if(nABs %% 20 > 0){nPlots<-nPlots+1}
for(i in 1:nPlots){
  range1=(i-1)*20+1
  range2=(i-1)*20+20
  pdf(paste0(figDir,"feature_plot_cite_harmony_",i,".pdf"),width=12,height=10)
  p<-FeaturePlot(citeCells, features = deABs[range1:range2], min.cutoff = "q05", max.cutoff = "q95", 
      ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.5)
  print(p)
  dev.off()
}

citeCells.subset <- subset(x = citeCells, idents = c("10","26","6","7","5","28","12","27"))
citeCells.subset <- RenameIdents(object = citeCells.subset,
                           "10" = "CD8_T",
                           "26" = "ILC_NK",
                           "6" = "CD4_T1",
                           "7" = "CD4_T2",
                           "5" = "B_Cell",
                           "28" = "Plasma_B",
                           "12" = "Myeloid",
                           "27" = "pDC")
citeCells<-citeCells.subset
sctranscormFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_umap_sctransform.Rdata")
if(!file.exists(sctranscormFile)){
  citeCells <- SCTransform(object = citeCells, verbose = TRUE)
  citeCells <- RunPCA(object = citeCells, verbose = TRUE,npcs = 100)
  citeCells <- JackStraw(object = citeCells, dims=100,num.replicate=100,verbose=T)
  citeCells <- ScoreJackStraw(citeCells,dims=1:100)
  save(citeCells,file=sctranscormFile)
} else {load(sctranscormFile)}



umapFile=paste0(homeDir,"project_results/Robjects/",sampleName,"_",treatment,"_umap_harmony.Rdata")
if(!file.exists(tSNEFile)){
  pcsData<-as.data.frame(citeCells@reductions$pca@jackstraw@overall.p.values)
  pcs<-pcsData[pcsData$Score<0.05,"PC"]
  citeCells <- FindNeighbors(object = citeCells, dims = pcs, verbose = FALSE)
  citeCells <- FindClusters(citeCells, dims = pcs, resolution = 0.8)
  citeCells <- RunUMAP(citeCells, dims = pcs)
  save(citeCells,file=umapFile)
} else {load(umapFile)}



adt.markers <- FindAllMarkers(citeCells, assay = "ADT", only.pos = TRUE)
deABs<-unique(adt.markers$gene)
deABs<-sort(deABs)
pdf(paste0(figDir,"CITE_clusters_merged_harmony2.pdf"),width=12,height=20)
p<-DoHeatmap(citeCells, features = unique(adt.markers$gene),raster=F, assay = "ADT", angle = 90) + NoLegend()
print(p)
dev.off()
  

pdf(paste0(figDir,"CITE_clusters_merged_harmony_tsne2.pdf"),width=12,height=8)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F","lightgreen","darkblue","lightblue","darkblue","orange","purple","gray","green","black","pink")
DimPlot(object = citeCells, reduction = "umap",
         label = TRUE, 
         pt.size = 2, 
         label.size = 6,
         cols = cols
)
dev.off()


genes<-c("nFeature_RNA","percent.mt",paste0(c("ACTB","CD14","TYR","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,"CITE_clusters_merged_harmony_feature_vst2.pdf"),width=12,height=16)
FeaturePlot(object = citeCells, features = genes, ncol = 3, pt.size = 0.5,cols = c("grey", "blue"), 
            reduction = "umap")
dev.off()

adt.markers <- FindAllMarkers(citeCells, assay = "ADT", only.pos = TRUE)
adtL<-split(adt.markers,adt.markers$gene)
minVal<-lapply(adtL,function(x){min(x[,1])})
minVal<-unlist(minVal)

deABs<-names(minVal)[order(minVal)]


#plot all cite data
nABs<-length(deABs)
nPlots<- nABs %/% 20
if(nABs %% 20 > 0){nPlots<-nPlots+1}
for(i in 1:nPlots){
  range1=(i-1)*20+1
  range2=(i-1)*20+20
  pdf(paste0(figDir,"feature_plot_cite_harmony2_",i,".pdf"),width=12,height=10)
  p<-FeaturePlot(citeCells, features = paste0("adt_",deABs[range1:range2]), min.cutoff = "q05", max.cutoff = "q95", 
      ncol = 4, cols = c("lightgrey", "blue"), reduction = "umap", pt.size = 0.2)
  print(p)
  dev.off()
}

df<-data.frame(ABs=names(sort(minVal)),p_value=as.numeric(sort(minVal)))
write.table(df,"MIA_abs_total.txt",quote=F,row.names=F,sep="\t")

write.table(adt.markers,"MIA_abs2.txt",quote=F,row.names=F,sep="\t")
adtL<-split(adt.markers,adt.markers$gene)



