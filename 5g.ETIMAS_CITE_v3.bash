#!/bin/bash

version=22
type="CITE"
sampleName="ETIMAS"$version
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sampleName"

#samples=( "3921_3963_4066" )
inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sampleName/"$type"_fastqs/"
whitelist=$inDir$sampleName"_whitelist.txt"
cellNum=`wc -l $whitelist | sed 's/ .*//'`
qsub -q short.q -hold_jid gzip -b y -j y -N $type$version -wd $logDir -pe smp 20 -V "CITE-seq-Count -T 20 -R1 "$inDir$version"_Spleen_"$type"_R1_001_combined.fastq.gz -R2 "$inDir$version"_Spleen_"$type"_R2_001_combined.fastq.gz -t "$inDir"tags.csv -cbf 1 -cbl 16 -umif 17 -umil 26 -o "$inDir$sampleName"_"$type".v3.tsv -cells "$cellNum" -wl "$whitelist
