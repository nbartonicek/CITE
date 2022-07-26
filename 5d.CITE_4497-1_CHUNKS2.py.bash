#!/bin/bash
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"

samples=( "4497-1_CHUNKS2" "4497-2_CELLSUS" )
#samples=( "4497-2_CELLSUSß" )

for sample in ${samples[@]};do
	rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sample"

	inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sample/CITE_fastqs/"
	cat $rawDir"/CITE_fastqs/raw/"*R1* > $inDir$sample"_CITESEQ_S1_R1_L001.fastq.gz"
	cat $rawDir"/CITE_fastqs/raw/"*R2* > $inDir$sample"_CITESEQ_S1_R2_L001.fastq.gz"

	#qsub -q short.q -hold_jid gzip -b y -j y -N CITE_$sample -wd $logDir -pe smp 12 -V "python /share/ScratchGeneral/nenbar/local/lib/CITE-seq-Count/CITE-seq-count.py -R1 $inDir$sample"_CITESEQ_S1_R1_L001.fastq.gz" -R2 $inDir$sample"_CITESEQ_S1_R2_L001.fastq.gz" -t $inDir"tags.csv" -cbf 1 -cbl 16 -umif 17 -umil 26 -tr \"^[ATGC]{15}[TGC][A]{6,}\" -hd 1 -o $inDir$sample"_CITE.tsv" -wl $inDir$sample\"_whitelist.txt\""
	qsub -q short.q -hold_jid gzip -b y -j y -N CITE_$sample -wd $logDir -pe smp 12 -V "python /share/ScratchGeneral/nenbar/local/lib/CITE-seq-Count/CITE-seq-count.py -R1 $inDir$sample"_CITESEQ_S1_R1_L001.fastq.gz" -R2 $inDir$sample"_CITESEQ_S1_R2_L001.fastq.gz" -t $inDir"tags.csv" -cbf 1 -cbl 16 -umif 17 -umil 26 -tr \"^[ATGC]{15}[TGC][A]{4,}\" -hd 1 -o $inDir$sample"_emptydrop_CITE.tsv" -wl $inDir$sample\"_emptydrop_whitelist.txt\""

done;
