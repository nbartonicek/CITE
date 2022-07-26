#!/bin/bash

module load gi/gcc/4.8.2
#raw files directory
homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/CITE"

scriptsPath="$homedir/projects/CITE/scripts"
logDir=$projectDir"/scripts/logs"
mkdir -p $logDir

orders=( "1" "3" "4" "5" )
orders=( "1" )
for order in ${orders[@]};do
	qsub -q short.q -b y -j y -N "c"$order -wd $logDir -pe smp 60 -P DSGClinicalGenomics -V "/share/ClusterShare/software/contrib/nenbar/src/R-3.5.2/builddir/bin/R --vanilla <$scriptsPath/2l.impute_Sunny.R --order=$order"
done;
