#!/bin/bash
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
inDir="/share/ScratchGeneral/nenbar/projects/CITE/sc_split/"
outDir="/share/ScratchGeneral/nenbar/projects/CITE/sc_split/vcfs/"

genomeFile="/share/ClusterShare/biodata/contrib/nenbar/genomes/human/cellranger/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa"
files=`ls /share/ScratchGeneral/nenbar/projects/CITE/sc_split/*REF_*.bam`

samtools view -S -b -q 10 -F 3844 3921_3963_4066_whitelist.bam > 3921_3963_4066_whitelist_filtered.bam
samtools rmdup 3921_3963_4066_whitelist_filtered.bam 3921_3963_4066_whitelist_filtered.rmdup.bam &


for file in ${files[@]};do
	echo $file
	sampleName=`basename $file | sed s/.*REF_//`
	outFile=`basename $file | sed s/.bam/.vcf/`
	qsub -q short.q -hold_jid gzip -b y -j y -N freebayes -wd $logDir -pe smp 5 -P DSGClinicalGenomics -V "freebayes -f $genomeFile -iXu -C 2 -q 1 $file >$outDir$outFile"
done;

inDir="/share/ScratchGeneral/nenbar/projects/CITE/sc_split/"
vcfFile=$inDir"all.vcf"
bamFile=$inDir"3921_3963_4066_whitelist_filtered.rmdup.sorted.bam"
barcodes=$inDir"3921_3963_4066_whitelist.txt"
refCount=$inDir"3921_3963_4066_refcount.txt"
altCount=$inDir"3921_3963_4066_altcount.txt"

qsub -q short.q -hold_jid gzip -b y -j y -N sc_split_matrices -wd $logDir -pe smp 20 -P DSGClinicalGenomics -V python $inDir"matrices.py" -v $vcfFile -i $bamFile -b $barcodes -r $refCount -a $altCount 

qsub -q short.q -hold_jid gzip -b y -j y -N sc_split_main -wd $logDir -pe smp 20 -P DSGClinicalGenomics -V python $inDir"main.py" -n 3 -r $refCount -a $altCount 

qsub -q short.q -hold_jid gzip -b y -j y -N sc_split_main -wd $logDir -pe smp 20 -P DSGClinicalGenomics -V python $inDir"genotype.py" -r $refCount -a $altCount 
