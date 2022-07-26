#!/bin/bash

type="HTO"
sampleNames=( "MIA" "pool1" "TONS_PBMC" )
sampleNames=( "MIA" )

for sampleName in ${sampleNames[@]};do
	logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
	rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sampleName"

#samples=( "3921_3963_4066" )
	inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sampleName/"$type"_fastqs/"
	whitelist=$inDir$sampleName"_whitelist.txt"
	cellNum=`wc -l $whitelist | sed 's/ .*//'`
	qsub -q short.q -hold_jid gzip -b y -j y -N $type$sampleName -wd $logDir -pe smp 12 -P DSGClinicalGenomics -V "python3 /share/ScratchGeneral/nenbar/local/lib/CITE-seq-Count/CITE-seq-count.py -R1 "$inDir$sampleName"_"$type"_R1_001.fastq.gz -R2 "$inDir$sampleName"_"$type"_R2_001.fastq.gz -t "$inDir"tags.csv -cbf 1 -cbl 16 -umif 17 -umil 28 -tr \"^[ATGC]{15}[TGC][ATGC]{4,}\" -hd 1 -o "$inDir$sampleName"_"$type".v3.tsv -wl "$whitelist 
done;
