#!/bin/bash
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"

samples=( "4290B" )
rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$samples"
inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$samples/HTO_fastqs/"
qsub -q short.q -hold_jid gzip -b y -j y -N HTO -wd $logDir -pe smp 12 -P DSGClinicalGenomics -V "python /share/ScratchGeneral/nenbar/local/lib/CITE-seq-Count/CITE-seq-count.py -R1 $inDir$samples"_HTO_S1_R1_001.fastq.gz" -R2 $inDir$samples"_HTO_S1_R2_001.fastq.gz" -t $inDir"all_tags.csv" -cbf 1 -cbl 16 -umif 17 -umil 26 -tr \"^[ATGC]{15}[TGCA][ATGC]{4,}\" -hd 2 -o $inDir$samples"_all.tsv" -wl $inDir$samples\"_whitelist.txt\""
