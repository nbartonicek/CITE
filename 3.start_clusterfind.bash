#!/bin/bash

module load briglo/R/3.4.2 
#module load gi/bcftools/1.0
#directory hierarchy
#raw files directory
homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/CITE"

scriptsPath="$homedir/projects/CITE/scripts"
logDir=$projectDir"/scripts/logs"
mkdir -p $logDir

dims=( "15" "20" "25" "30" "35" "40" )
dims=( "40" "50" )
resolutions=( "0.6" "0.8" )
#resolutions=( "1" "1.2" )
perplexities=( "5" "10" "15" "20" "50")

#dims=( "15" )
#resolutions=( "0.6" )
#perplexities=( "5" )
for dim in ${dims[@]};do
	for resolution in ${resolutions[@]};do
		for perplexity in ${perplexities[@]};do
			echo $dim $resolution $perplexity
			qsub -q long.q -b y -j y -N "c"$dim.$resolution.$perplexity -wd $logDir -pe smp 5 -V "R --vanilla <$scriptsPath/3.find_clusters.R --dims=$dim --resolution=$resolution --perplexity=$perplexity"
		done;

	done;
done;
#data=read.table("test.tsv",header=T,row.names=1,sep=",")
#apply(data,1,function(x){sum(x>1)})
