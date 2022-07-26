#!/bin/bash
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4530"

samples=( "CID4530_CITE" "CID4530_HTO" "CID4530-N_CITE" "CID4530-N_HTO" )
inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4515/CITE_fastqs/"
qsub -q short.q -hold_jid gzip -b y -j y -N CITE -wd $logDir -pe smp 12 -P DSGClinicalGenomics -V "python /share/ScratchGeneral/nenbar/local/lib/CITE-seq-Count/CITE-seq-count.py -R1 $inDir"CID4515_CITESEQ_S1_R1_L001.fastq.gz" -R2 $inDir"CID4515_CITESEQ_S1_R2_L001.fastq.gz" -t $inDir"tags.csv" -cbf 1 -cbl 16 -umif 17 -umil 26 -tr \"^[ATGC]{15}[TGC][A]{6,}\" -hd 1 -o $inDir"4515_CITE.tsv" -wl $inDir\"whitelist.txt\""
