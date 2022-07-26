#!/bin/bash

type="CITE"
sampleNames=( "4497-2" "4515" "SVH130918A_control" "3921_3963_4066" "4497-1" "4290B" "4530-N" "4530-T" "4T1_3x" "ETIMAS_S21" "ETIMAS_S22" "ETIMAS_T21" "ETIMAS_T22" "MIA" "pool1" "TONS_PBMC" )
sampleNames=( "evaapo_02" )
method="umis"
for sampleName in ${sampleNames[@]};do
	logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
	rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sampleName"

#samples=( "3921_3963_4066" )
	inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sampleName/"$type"_fastqs/"
	whitelist=$inDir$sampleName"_whitelist.txt"
	cellNum=`wc -l $whitelist | sed 's/ .*//'`

	seqlengths=( "81" "79" "77" "75" "73" )
	#for seqlength in ${seqlengths[@]};do
	R1=$inDir$sampleName"_"$type"_R1_001.fastq.gz"
	R2=$inDir$sampleName"_"$type"_R2_001.fastq.gz"
	R1Length=`zless $R1 | sed 1d | head -n1 | wc -m `
	echo $R1
	echo $R1Length
	umil=28
	barcode="CCCCCCCCCCCCCCCCNNNNNNNNNNNN"


	if [ $R1Length == "27" ]; then
		umil="26"
		barcode="CCCCCCCCCCCCCCCCNNNNNNNNNN"
	fi

	qsub -q short.q -hold_jid gzip -b y -j y -N $type$sampleName -wd $logDir -pe smp 50 -V "CITE-seq-Count -T 25 -R1 "$R1" -R2 "$R2" -t "$inDir"CITE_ListV3_157ADT.csv -cbf 1 -cbl 16 -umif 17 -umil $umil -o "$inDir$sampleName"_"$type".v3 -cells "$cellNum" -wl "$whitelist 
	method="reads"

	#qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp 4 -P DSGClinicalGenomics -V"
	#$qsubLine -N whitelist$sampleName "umi_tools whitelist --stdin "$R1" \
	#	--plot-prefix=$inDir/umiTools_"$method" --set-cell-number="$cellNum" \
	#	--extract-method=string --method="$method" --bc-pattern="$barcode"  \
	#	--log2stderr > $inDir/umitools_bcs_reads.txt;"
	#method="umis"
	#$qsubLine -N whitelist$sampleName "umi_tools whitelist --stdin "$R1" \
	#		--plot-prefix=$inDir/umiTools_"$method" --set-cell-number="$cellNum" \
	#		--extract-method=string --method="$method" --bc-pattern="$barcode"  \
	#		--log2stderr > $inDir/umitools_bcs_UMIs.txt;"

	#done;
done;
