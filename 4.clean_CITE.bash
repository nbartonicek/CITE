#!/bin/bash
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
qsub -q short.q -b y -j y -N gzip -wd $logDir -pe smp 8 -V "zless /share/ScratchGeneral/nenbar/projects/CITE/raw_files/4515/CITE_fastqs/CID4515_S2_R2_001.fastq.gz |  sed 's/.*CACCCGAGAATTCCA//' | gzip -c > /share/ScratchGeneral/nenbar/projects/CITE/raw_files/4515/CITE_fastqs/clean_R2.fastq.gz"