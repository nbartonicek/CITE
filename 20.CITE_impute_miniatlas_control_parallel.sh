#!/bin/bash

projectname="CITE"
sampleNames=( "CID3838"  "CID4040"  "CID4515"  "CID4378N" "CID3946" )
sampleNames=( "CID4515" )
scriptsPath="/share/ScratchGeneral/nenbar/projects/$projectname/scripts"
logDir="$scriptsPath/logs"
mkdir -p $logDir

numcores=40

for sampleName in ${sampleNames[@]}; do
	qsub -q short.q -hold_jid gzip -b y -j y -N CITE_$sampleName -wd $logDir -pe smp $numcores -P TumourProgression -V "R --vanilla <$scriptsPath/20.CITE_impute_miniatlas_control.R --eliminateSample=$sampleName"
done;