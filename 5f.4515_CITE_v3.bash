#!/bin/bash
sampleName="4515"
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sampleName"

#samples=( "3921_3963_4066" )
inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sampleName/CITE_fastqs/"
qsub -q short.q -hold_jid gzip -b y -j y -N CITE -wd $logDir -pe smp 20 -P DSGClinicalGenomics -V "CITE-seq-Count -T 20 -R1 "$inDir"CID"$sampleName"_CITESEQ_S1_R1_L001.fastq.gz" -R2 "$inDir"CID"$sampleName"_CITESEQ_S1_R2_L001.fastq.gz" -t $inDir"tags.csv" -cbf 1 -cbl 16 -umif 17 -umil 26 -o $inDir$sampleName"_CITE.v3.tsv" -cells 6633 -wl $inDir$sampleName\"_whitelist.txt\""
