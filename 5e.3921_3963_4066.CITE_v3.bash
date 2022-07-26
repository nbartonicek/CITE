#!/bin/bash
sampleName="3921_3963_4066"
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sampleName"

inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sampleName/CITE_fastqs/"

cellNum=`wc -l $inDir$sampleName"_whitelist.txt" | sed 's/ .*//'`
qsub -q short.q -hold_jid gzip -b y -j y -N CITE -wd $logDir -pe smp 60 -P DSGClinicalGenomics -V "CITE-seq-Count -T 60 -R1 "$inDir$sampleName"_CITE_S4_R1_001.fastq.gz" -R2 $inDir$sampleName"_CITE_S4_R2_001.fastq.gz" -t $inDir"tags.csv -cbf 1 -cbl 16 -umif 17 -umil 26 -o "$inDir$sampleName"_CITE.v3.tsv -cells $cellNum -wl $inDir$sampleName\"_whitelist.txt\""
