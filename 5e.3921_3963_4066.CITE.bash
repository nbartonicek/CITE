#!/bin/bash
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/3921_3963_4066"

samples=( "3921_3963_4066" )
inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/3921_3963_4066/CITE_fastqs/"
qsub -q short.q -hold_jid gzip -b y -j y -N CITE -wd $logDir -pe smp 12 -P DSGClinicalGenomics -V "python /share/ScratchGeneral/nenbar/local/lib/CITE-seq-Count/CITE-seq-count.py -R1 $inDir"3921_3963_4066_CITE_S4_R1_001.fastq.gz" -R2 $inDir"3921_3963_4066_CITE_S4_R2_001.fastq.gz" -t $inDir"all_tags.csv" -cbf 1 -cbl 16 -umif 17 -umil 26 -tr \"^[ATGC]{15}[TGC][ATGC]{4,}\" -hd 1 -o $inDir"3921_3963_4066_CITE.tsv" -wl $inDir\"3921_3963_4066_whitelist.txt\""
