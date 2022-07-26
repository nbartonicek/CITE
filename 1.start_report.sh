#!/bin/bash
module load gi/gcc/4.8.2
module load briglo/R/3.4.2 
#export PIPELINE_PATH="/share/ScratchGeneral/nenbar/projects/single_cell"
#The analysis is based on UMItools: https://github.com/CGATOxford/UMI-tools/blob/master/doc/Single_cell_tutorial.md
#It consists of a set of bash job submissions of individial methods

#Input files: 
#file: gtf of annotation
#file: fa of genome
#file: fa of transcriptome
#file: flat transcriptomes that prevent false multimaper loss (cgat) https://github.com/CGATOxford/cgat
#directory: star index based on flat transcriptomes (if available, above 4 are redundant)
#string: sample name
#string: platform
#string: protocol version
#string: organism
#directory: output files
#directory: logs
#string: group on cluster
#integer: number of cores for a small job
#integer: number of cores for a large job
#string: analysis parameters (add)


#### preparing additional files #####
#### prepare transcriptomes
#take gencode annotation, and keep only annotations with "protein_coding" "antisense_RNA"      "lincRNA"
#data=import("gencode.vM16.annotation.gtf")
##watch out, it is MT for human and mt for mouse
#mitoGenes<-unique(data$gene_name[grepl("^MT-",data$gene_name)])
#dataM<-data[data$gene_name %in% mitoGenes]
##since mito genes are spread across gene_type category,have to first remove, then add them
#dataT=data[data$gene_type %in% c("protein_coding","antisense_RNA","lincRNA")]
#dataT=dataT[!(dataT$gene_name %in% mitoGenes)]
#dataT=c(dataT,dataM)
#export.gff2(dataT,"gencode.v27.polyA.annotation.gtf")
###to make the extended version do this

#library(rtracklayer)
#data=import("gencode.vM16.polyA.annotation.gtf")
#data$id=1:length(data)
#dataS<-data
#values(dataS)<-NULL
#strands<-strand(data)
#strand(dataS)="+"
#dataS$id=data$id
#dataL<-split(dataS,data$gene_id)
##dataLL<-endoapply(dataL[1:1000],function(x){c(x[1],x[length(x)]);cat(".")})
#dataLS<-endoapply(dataL,function(x){x[1]=resize(x[1],width(x[1])+500,fix="end");cat(".");num=length(x);x[num]=resize(x[num],width(x[num])+500);x})
#dataUL<-unlist(dataLS)
#dataUL<-dataUL[order(dataUL$id)]
#strand(dataUL)<-strands
#export.gff(dataUL,"gencode.vM16.polyA.extended.annotation.gtf")
#source /share/ScratchGeneral/nenbar/local/lib/CGAT/conda-install/bin/activate cgat-s
#cgat gtf2gtf --method=merge-exons -I gencode.v27.polyA.annotation.gtf | cgat gtf2gtf --method=set-transcript-to-gene -S genes_merged.gtf
#source deactivate
###### fasta file
#module load gi/boost/1.53.0
#module load gi/cufflinks/2.2.1
#gffread -w genes_merged.fa -g ../fasta/genome.fa genes_merged.gtf
#gffread -w genes_merged.fa -g /share/ClusterShare/biodata/contrib/nenbar/genomes/star/mm10/mm10.fa genes_merged.gtf
#gffread -w genes_merged.fa -g /share/ClusterShare/biodata/contrib/nenbar/genomes/hg38_mm10/hg38_mm10.fa genes_merged.gtf
#### build star index
#/share/ClusterShare/biodata/contrib/nenbar/projects/single_cell/pipelines/cellranger-2.0.2/STAR/2.5.1b/STAR  --runMode genomeGenerate --genomeDir "." --genomeFastaFiles "genes_merged.fa" --genomeChrBinNbits 12

#this is loading human and mouse index, but that should be taken from the manifest file in the future



############ information for running jobs on the cluster ##########
scriptDir="$PIPELINE_PATH/scripts/"
#genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/mm10_gencode.v16"
logDir=$scriptDir"logs"
numcoresSmall=10
queue="short"
mkdir -p $logDir
tag="-P DSGClinicalGenomics"
#tag=""

############ the only required piece of information when the correct samplesheet available in projects dir ##########

projectnames=( "201802_DAVGAL_10X_01" "180310_uyengu_Indrop_Davgal_pool3"  "20180301_uyengu_Dropseq_Davgal" "20180226_uyengu_Indrop_Davgal" )
projectnames=( "201802_DAVGAL_10X_01" "180310_uyengu_Indrop_Davgal_pool3" "20180301_uyengu_Dropseq_Davgal" )
#basespaceID="71326261"
for projectname in ${projectnames[@]};do
	echo $projectname

	filename=$scriptDir"/projects/"$projectname".csv"

	unset cellnumbers
	unset speciesList
	declare -A cellnumbers
	declare -A speciesList
	declare -A replicateList
	{
		read;
		while read -r -a line || [ -n "$line" ]; do
			#echo $line
	    	platform=`echo $line | cut -f1 -d ","`
	    	sample_name=`echo $line | cut -f3 -d ","`
	    	species=`echo $line | cut -f4 -d ","`
			cellnumber=`echo $line | cut -f6 -d ","`
			replicate_name=`echo $line | cut -f9 -d ","`
			echo $replicate_name

	    	cellnumbers[$replicate_name]=$((${cellnumbers[$replicate_name]}+$cellnumber));
			speciesList[$replicate_name]=$species
		done 
	}< "$filename"


	for sample in ${!cellnumbers[@]}; do
		echo ${cellnumbers["$sample"]}
		#run the script 6.make_reports.R which makes R objects with data for reports
		qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp $numcoresSmall $tag -V"
		$qsubLine -N rpt$projectname "Rscript "$scriptDir"6.make_reports.R "$projectname" "$sample" "$platform" "${speciesList["$sample"]}" "${cellnumbers["$sample"]}" T"
		#line=`echo $qsubLine -N rpt$projectname "Rscript "$scriptDir"6.make_reports.R "$projectname" "$sample" "$platform" "${speciesList["$sample"]}" "${cellnumbers["$sample"]}" T"`
		#echo $line


	done;
done;
