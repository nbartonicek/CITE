#!/bin/bash
module load gi/seqtk/1.0
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530"

samples=( "CID4530_" "CID4530-N_" )
#samples=( "CID4530-N_" )
inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530/CITE_fastqs"
outDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530/CITE_fastqs/subsamples"
mkdir -p $inDir
mkdir -p $outDir

classes=( "CITE" "HTO" )
#classes=( "CITE" )
reads=( "R1" "R2" )

subsamples={1..9}

decimal="00"

for sample in ${samples[@]}; do
	for class in ${classes[@]}; do
		echo $sample$class
		for replicate in {1..3}; do
			for subsample in {1..9}; do
				for read in ${reads[@]}; do
					qsub -q short.q -b y -j y -N seqtk$decimal$sample$subsample -wd $logDir -pe smp 2 -P DSGClinicalGenomics -V "seqtk sample -s10$replicate $inDir/$sample$class"_"$read".fastq.gz" 0.$decimal$subsample | gzip -c > $outDir/$sample$class"_"$read"_0.$decimal$subsample.$replicate.fastq.gz""
				done;
				qsub -q short.q -b y -j y -hold_jid seqtk$decimal$sample$subsample -N CITE$decimal$sample$subsample -wd $logDir -pe smp 6 -P DSGClinicalGenomics -V "python /share/ScratchGeneral/nenbar/local/lib/CITE-seq-Count/CITE-seq-count.py -R1 $outDir/$sample$class"_R1_0.$decimal$subsample.$replicate.fastq.gz" -R2 $outDir/$sample$class"_R2_0.$decimal$subsample.$replicate.fastq.gz" -t $inDir/$class"_tags.csv" -cbf 1 -cbl 16 -umif 17 -umil 26 -tr \"^[ATGC]{15}[TGC][A]{6,}\" -hd 1 -o $outDir/$sample$class"_0.$decimal$subsample.$replicate.emptydrop.tsv" -wl $inDir/$sample\"emptydrop_whitelist.txt\""
			done;
		done;
	done;
done;
