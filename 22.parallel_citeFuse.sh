#!/bin/bash

projectname="CITE"
inDir="/share/ScratchGeneral/nenbar/projects/CITE_paper/project_results/Robjects/citefuse"

files=`ls $inDir | grep 'umap'`


#sampleNames=( "4378N" )
scriptsPath="/share/ScratchGeneral/nenbar/projects/$projectname/scripts"
logDir="$scriptsPath/logs"
mkdir -p $logDir

numcores=20

for file in ${files[@]}; do
	qsub -q short.q -b y -j y -N CITE_$file -wd $logDir -pe smp $numcores -P TumourProgression -V "R --vanilla <$scriptsPath/22.citefuse.test.R --fileName=$file"
done;