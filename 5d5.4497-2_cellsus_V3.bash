#!/bin/bash
sampleName="4497-2_CELLSUS"
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sampleName"

inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sampleName/CITE_fastqs/"

cellNum=`wc -l $inDir$sampleName"_whitelist.txt" | sed 's/ .*//'`
numcor=40
qsub -q short.q -hold_jid gzip -b y -j y -N CITE -wd $logDir -pe smp $numcor -P DSGClinicalGenomics -V "CITE-seq-Count -T $numcor -R1 "$inDir$sampleName"_CITESEQ_S1_R1_L001.fastq.gz" -R2 $inDir$sampleName"_CITESEQ_S1_R2_L001.fastq.gz" -t $inDir"tags.csv -cbf 1 -cbl 16 -umif 17 -umil 26 -o "$inDir$sampleName"_CITE.v3.tsv -cells $cellNum -wl $inDir$sampleName\"_whitelist.txt\""
