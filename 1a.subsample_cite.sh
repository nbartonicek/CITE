#!/bin/bash
module load gi/seqtk/1.0
inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files"
projectname="4378N"
logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
#mkdir -p $inDir/$projectname\.10/CITE_fastqs
#cp $inDir/$projectname/CITE_fastqs/*.gz $inDir/$projectname\.10/CITE_fastqs
qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp 5 -V"	
	$qsubLine -N gz gunzip $inDir/$projectname\.10/CITE_fastqs/*R1*.gz
	$qsubLine -N gz gunzip $inDir/$projectname\.10/CITE_fastqs/*R2*.gz
for i in {1..9};do 
	mkdir -p $inDir/$projectname\.$i/CITE_fastqs
	#qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp 5 -P TumourProgression -V"	
	$qsubLine -hold_jid gz -N CITE$i.R1 "seqtk sample -s100 $inDir/$projectname\.10/CITE_fastqs/$projectname\_CITE_R1_001.fastq 0.$i | gzip -c >$inDir/$projectname\.$i/CITE_fastqs/$projectname\_CITE_R1_001.fastq.gz"
	$qsubLine -hold_jid gz -N CITE$i.R2 "seqtk sample -s100 $inDir/$projectname\.10/CITE_fastqs/$projectname\_CITE_R2_001.fastq 0.$i | gzip -c >$inDir/$projectname\.$i/CITE_fastqs/$projectname\_CITE_R2_001.fastq.gz"
done 
