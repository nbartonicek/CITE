#!/bin/bash

type="CITE"
sampleNames=( "3838"  "4040"  "4515"  "3946"  "4378N" )
method="umis"
for sampleName in ${sampleNames[@]};do
	logDir="/share/ScratchGeneral/nenbar/projects/CITE/scripts/logs"
	rawDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/CITE_fastqs/$sampleName"
	tagDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/"
#samples=( "3921_3963_4066" )
	inDir="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/$sampleName/CITE_fastqs/"
	whitelist=$inDir$sampleName"_whitelist_miniatlas.txt"
	
	#gunzip -c /share/ScratchGeneral/sunwu/cellranger/cellranger_count_files_temp/GA_mRNA+LB_*-1*092019/run01_sunwupipeline/output/count/CID$sampleName/count_CID"$sampleName"_GRCh38_mm10/outs/raw_feature_bc_matrix/barcodes.tsv.gz > $whitelist

	cellNum=`wc -l $whitelist | sed 's/ .*//'`


	R1=$inDir$sampleName"-"$type"_R1_001.fastq.gz"
	R2=$inDir$sampleName"-"$type"_R2_001.fastq.gz"
	R1Length=`zless $R1 | sed 1d | head -n1 | wc -m `

	echo $R1Length
	umil=28
	barcode="CCCCCCCCCCCCCCCCNNNNNNNNNNNN"

	#echo $R1Length
	#umil=24
	#barcode="CCCCCCCCCCCCCCCCNNNNNNNN"

	if [ $R1Length == "27" ]; then
		umil="26"
		barcode="CCCCCCCCCCCCCCCCNNNNNNNNNN"
	fi

	qsub -q short.q -hold_jid gzip -b y -j y -N $type$sampleName -wd $logDir -pe smp 30 -P TumourProgression -V "CITE-seq-Count -T 10 -R1 "$R1" -R2 "$R2" -t "$tagDir"tags_118.csv -cbf 1 -cbl 16 -umif 17 -umil $umil -o "$inDir$sampleName"_"$type".miniatlas -cells "$cellNum" -wl "$whitelist 
	#method="reads"

	#qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp 4 -P TumourProgression -V"
	#$qsubLine -N whitelist$sampleName$method "umi_tools whitelist --stdin "$R1" \
	#	--plot-prefix=$inDir/umiTools_"$method" --set-cell-number=5000 \
	#	--extract-method=string --method="$method" --bc-pattern="$barcode"  \
	#	--log2stderr > $inDir/umitools_bcs_reads.txt;"
	#method="umis"
	#$qsubLine -N whitelist$sampleName$method "umi_tools whitelist --stdin "$R1" \
	#		--plot-prefix=$inDir/umiTools_"$method" --set-cell-number=5000 \
	#		--extract-method=string --method="$method" --bc-pattern="$barcode"  \
	#		--log2stderr > $inDir/umitools_bcs_UMIs.txt;"
#
	#qsub -q short.q -hold_jid whitelist$sampleName$method -b y -j y -N $type$sampleName -wd $logDir -pe smp 20 -V "cut -f1 $inDir/umitools_bcs_UMIs.txt > $whitelist; CITE-seq-Count -T 10 -R1 "$R1" -R2 "$R2" -t "$inDir"tags.csv -cbf 1 -cbl 16 -umif 17 -umil $umil -o "$inDir$sampleName"_"$type".v3 -cells 8361 -wl "$whitelist 
	#method="reads"


done;
