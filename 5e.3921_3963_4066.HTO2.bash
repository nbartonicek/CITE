#!/bin/bash
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/3921_3963_4066"

samples=( "3921_3963_4066" )
inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/3921_3963_4066/HTO_fastqs/"
qsub -q short.q -hold_jid gzip -b y -j y -N HTO -wd $logDir -pe smp 40 -P DSGClinicalGenomics -V "CITE-seq-Count -T 40 -R1 $inDir"3921_3963_4066_HTO_S1_R1_001.fastq.gz" -R2 $inDir"3921_3963_4066_HTO_S1_R2_001.fastq.gz" -t $inDir"all_tags.csv" -cbf 1 -cbl 16 -umif 17 -umil 26 --max-errors 2 -o $inDir"3921_3963_4066_highMem.tsv" -cells 13101 -wl $inDir\"3921_3963_4066_whitelist.txt\""
