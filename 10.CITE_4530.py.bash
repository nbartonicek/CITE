#!/bin/bash
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530"

samples=( "CID4530_" "CID4530-N_" )
#samples=( "CID4530-N_" )
inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530/CITE_fastqs"
mkdir -p $inDir
classes=( "CITE" "HTO" )
classes=( "CITE" )
reads=( "R1" "R2" )

for sample in ${samples[@]}; do
	for class in ${classes[@]}; do
		echo $sample$class
		#for read in ${reads[@]}; do
		#	cat $rawDir/$sample$class*$read* >$inDir/$sample$class"_"$read".fastq.gz"
		#done;
		qsub -q short.q -b y -j y -N CITE$sample -wd $logDir -pe smp 12 -P DSGClinicalGenomics -V "python /share/ScratchGeneral/nenbar/local/lib/CITE-seq-Count/CITE-seq-count.py -R1 $inDir/$sample$class"_R1.fastq.gz" -R2 $inDir/$sample$class"_R2.fastq.gz" -t $inDir/$class"_tags.csv" -cbf 1 -cbl 16 -umif 17 -umil 26 -tr \"^[ATGC]{15}[TGC][A]{6,}\" -hd 1 -o $inDir/$sample$class"_emptydrop.tsv" -wl $inDir/$sample\"emptydrop_whitelist.txt\""
	done;
done;
