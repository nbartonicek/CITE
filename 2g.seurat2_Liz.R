library(Seurat)
library(dplyr)
library(Matrix)
library(devtools)
library(sctransform)
library(ggplot2)

homeDir="../../../projects/Liz/"
figDir=paste0(homeDir,"project_results/figures/")
projectname="Liz"
sampleName="190121_alldat_ano_pseudo_cycleBasecc"
treatment="v2"
#treatment="mt0.5_nGene7000"
figDir=paste0(homeDir,"project_results/figures/",sampleName,"_",treatment,"/")
system(paste0("mkdir -p ",figDir))

tSNEFile=paste0(homeDir,"/raw_files/",sampleName,".rdata")
load(tSNEFile)


pdf(paste0(figDir,sampleName,"_",treatment,"_AA_vs_PA.pdf"),width=12,height=8)
cols<-c("#b2182b", "#1F618D", "#053061", 
                        "#bebada", "#f4a582", "#1c9099", 
                        "#85929E", "#9B59B6", "#74add1", 
                        "#de77ae", "#76D7C4", "#b8e186", 
                        "#4393c3", "#c51b7d", "#b8e186", 
                        "#17A589", "#A93226", "#F5B041", "#F1C40F","darkblue","orange")
DimPlot(object = alldat, dim1=Pseudotime,dim2=x_pseudo, 
         do.label = T, 
         pt.size = 2, 
         label.size = 6
)
dev.off()




genes<-c("nFeature_RNA","percent.mt",paste0(c("ACTB","CD14","EPCAM","PDGFRA","FCGR3A","COL1A1","CD19","CD3E","CD4","CD8A","KRT8","KRT5", "PECAM1","PTPRC", "CD34", "ERBB2")))
pdf(paste0(figDir,sampleName,"_",treatment,"_feature.pdf"),width=12,height=16)
FeaturePlot(object = alldat, features.plot(c(x_pseudo,y_pseudo)), features = genes, nCol = 3, pt.size = 0.5,cols = c("grey", "blue"))

dev.off()
#}



#add CITE data
cite <-cite[,colnames(cite) %in% colnames(citeCells@data) ]
citeCells <- SetAssayData(citeCells, assay.type = "CITE", slot = "raw.data", new.data = cite)
citeCells <- NormalizeData(citeCells, assay.type = "CITE", normalization.method = "genesCLR")
citeCells <- ScaleData(citeCells, assay.type = "CITE", display.progress = FALSE)

pdf(paste0(figDir,"feature_plot_cite_",sampleName,".pdf"),width=12,height=20)
FeaturePlot(citeCells, features.plot = row.names(cite), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 4, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


pdf(paste0(figDir,"ridgeplot_",sampleName,".pdf"),width=12,height=4)

RidgePlot(citeCells, features.plot = row.names(cite), cols.use = cols,nCol = 3)

dev.off()


genes<-c("CD14","CD16","COL1A1","EPCAM","CD19","CD3","","","")
#testing for dead cells
FeaturePlot(citeCells, features.plot = c("MT-ND1","mm10___mt-Nd1","MYC","MGAT1"), min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 3, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
save(citeCells,file="../project_results/Robjects/4497_PBMC_CITE.Rdata")


cbmc[["ADT"]] <- CreateAssayObject(counts = cite)
cbmc <- NormalizeData(cbmc, assay = "ADT", normalization.method = "CLR")
cbmc <- ScaleData(object = cbmc, assay = "ADT")
FeaturePlot(object = cbmc, features = c("adt_CD3", "adt_CD11c", "adt_CD8", "adt_CD16", "CD3E", "ITGAX", 
    "CD8A", "FCGR3A"), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)

