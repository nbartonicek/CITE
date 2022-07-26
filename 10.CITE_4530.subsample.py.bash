#!/bin/bash
module load gi/seqtk/1.0
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530"

samples=( "CID4530_" "CID4530-N_" )
samples=( "CID4530-N_" )
inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530/CITE_fastqs"
outDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530/CITE_fastqs/subsamples"
mkdir -p $inDir
mkdir -p $outDir

classes=( "CITE" "HTO" )
classes=( "CITE" )
reads=( "R1" "R2" )

subsamples={1..9}

for sample in ${samples[@]}; do
	for class in ${classes[@]}; do
		echo $sample$class
		for subsample in {1..2}; do
			seqtk sample -s100 $inDir/$sample$class"_"$read".fastq.gz" 0.$subsample | gzip -c > $outDir/$sample$class"_"$read"_0.$subsample.fastq.gz"
			qsub -q short.q -b y -j y -N CITE$sample$subsample -wd $logDir -pe smp 12 -P DSGClinicalGenomics -V "python /share/ScratchGeneral/nenbar/local/lib/CITE-seq-Count/CITE-seq-count.py -R1 $outDir/$sample$class"_R1_0.$subsample.fastq.gz" -R2 $outDir/$sample$class"_R2_0.$subsample.fastq.gz" -t $inDir/$class"_tags.csv" -cbf 1 -cbl 16 -umif 17 -umil 26 -tr \"^[ATGC]{15}[TGC][A]{6,}\" -hd 1 -o $outDir/$sample$class"_0.$subsample.tsv" -wl $inDir/$sample\"whitelist.txt\""
		done;
	done;
done;
