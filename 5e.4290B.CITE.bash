#!/bin/bash
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4290B"

samples=( "4290B" )
inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4290B/CITE_fastqs/"
qsub -q short.q -hold_jid gzip -b y -j y -N CITE -wd $logDir -pe smp 12 -P DSGClinicalGenomics -V "python /share/ScratchGeneral/nenbar/local/lib/CITE-seq-Count/CITE-seq-count.py -R1 $inDir"4290B_CITE_S1_R1_001.fastq.gz" -R2 $inDir"4290B_CITE_S1_R2_001.fastq.gz" -t $inDir"all_tags.csv" -cbf 1 -cbl 16 -umif 17 -umil 26 -tr \"^[ATGC]{15}[TGC][ATGC]{4,}\" -hd 1 -o $inDir"4290B_CITE.tsv" -wl $inDir\"4290B_whitelist.txt\""
