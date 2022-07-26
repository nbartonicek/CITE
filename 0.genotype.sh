#!/bin/bash

module load gi/gcc/4.8.2
module load gi/picard-tools/1.138
#raw files directory
homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/CITE"

scriptsPath="$homedir/projects/CITE/scripts"
logDir=$projectDir"/scripts/logs"
mkdir -p $logDir

projectname="3921_3963_4066_Blood"

annotationDir=$projectDir/"genotyping/annotation"
inDir=$projectDir"/raw_files/genotype/YAN6909_UKB_A01-F10_RESULTS/"
outDir=$projectDir/"genotyping/$projectname"
step1Dir=$outDir/step1
step2Dir=$outDir/step2

mkdir -p $outDir
mkdir -p $step1Dir
mkdir -p $step2Dir


(echo cel_files; ls -1 $inDir/*.CEL) > $step1Dir/"cel_list1.txt"

############### remove the samples that have failed - with a DQC value less than the default DQC threshold of 0.82
 
more $inDir"YAN6909_UKB_1095_BP_Workflow_QC_Report_Table_Arrays_A01-F10.csv" | grep "Fail" | cut -f1 -d "," | grep "CEL" >$step1Dir/"failed.txt"

(echo cel_files; grep $step1Dir/"cel_list1.txt" -v -f $step1Dir/"failed.txt" | grep "Blood" | grep '3921\|3963\|4066') >$step1Dir/"cel_list2.txt"

############### run apt-genotype-axiom

#/share/ScratchGeneral/nenbar/local/lib/apt-2.10.2.2-x86_64-intel-linux/bin/apt-genotype-axiom --log-file $inDir"apt-genotype-axiom.log" --arg-file $inDir"Axiom_UKB_WCSG.r5.apt-genotype-axiom.AxiomCN_GT1.apt2.xml" --analysis-files-path $inDir --out-dir $outDir --dual-channel-normalization true --cel-files $inDir/cel_list2.txt

genotypeLine="/share/ScratchGeneral/nenbar/local/lib/apt-2.10.2.2-x86_64-intel-linux/bin/apt-genotype-axiom --log-file $step1Dir/"apt-genotype-axiom.log" --arg-file $annotationDir/"Axiom_UKB_WCSG.r5.apt-genotype-axiom.AxiomCN_GT1.apt2.xml" --analysis-files-path $annotationDir --out-dir $step1Dir --dual-channel-normalization true --cel-files $step1Dir/"cel_list2.txt" --summaries --write-models --batch-folder $step1Dir"

qsub -q short.q -b y -j y -N genotypeStep1 -wd $logDir -pe smp 20 -P TumourProgression -V $genotypeLine


############### run the command line version of SNPolisher

psmetricsLine="/share/ScratchGeneral/nenbar/local/lib/apt-2.10.2.2-x86_64-intel-linux/bin/ps-metrics --posterior-file $step1Dir/AxiomGT1.snp-posteriors.txt --call-file $step1Dir/AxiomGT1.calls.txt --metrics-file $step1Dir/SNPolisher/metrics.txt" 

psclassLine="/share/ScratchGeneral/nenbar/local/lib/apt-2.10.2.2-x86_64-intel-linux/bin/ps-classification --species-type human --metrics-file $step1Dir/SNPolisher/metrics.txt --output-dir $step1Dir/SNPolisher --ps2snp-file $annotationDir/Axiom_UKB_WCSG.r5.ps2snp_map.ps"

otvLine="/share/ScratchGeneral/nenbar/local/lib/apt-2.10.2.2-x86_64-intel-linux/bin/otv-caller --pid-file $step1Dir/SNPolisher/Recommended.ps --batch-folder $step1Dir --output-dir $step1Dir/OTV"

qsub -q short.q -b y -j y -N genotypeStep2 -wd $logDir -pe smp 20 -P TumourProgression -V $psmetricsLine

qsub -q short.q -b y -j y -N genotypeStep3 -wd $logDir -pe smp 20 -P TumourProgression -V $psclassLine

qsub -q short.q -b y -j y -N genotypeStep4 -wd $logDir -pe smp 20 -P TumourProgression -V $otvLine


############### fetch sqlite annotation from https://www.thermofisher.com/order/catalog/product/901153?SID=srch-srp-901153 

outFolder="/share/ScratchGeneral/nenbar/projects/CITE/genotyping/$projectname/step2"

/share/ScratchGeneral/nenbar/local/lib/apt-2.10.2.2-x86_64-intel-linux/bin/apt-format-result --batch-folder $step1Dir/OTV --snp-list-file $step1Dir/OTV/OTV.keep.ps --annotation-file $annotationDir/Axiom_UKB_WCSG.na35.annot.db --export-vcf-file $outFolder"/"$projectname".b37.vcf"
java -jar /share/ClusterShare/software/contrib/gi/picard-tools/1.138/picard.jar CreateSequenceDictionary R=/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg38/hg38.fa O=/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg38/hg38.dict
more $outFolder"/"$projectname".b37.vcf" | grep -v UNKNOWNPOSITION >$outFolder"/"$projectname".b37.clean.vcf"
java -jar ~/local/bin/picard.jar LiftoverVcf I=$outFolder"/"$projectname".b37.clean.vcf" O=$outFolder"/"$projectname".hg38.vcf" CHAIN=$annotationDir"/b37ToHg38.over.chain" REJECT=$annotationDir"/rejected_variants.vcf" R=/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg38/hg38.fa

awk '{gsub(/^chr/,""); print}' $outFolder"/"$projectname".hg38.vcf" > $outFolder"/"$projectname".b38.vcf"
############### with or without splink

qsub -q short.q -b y -j y -N genotypeStep4 -wd $logDir -pe smp 20 -P TumourProgression -V demuxlet --alpha 0.5 --sam /share/ScratchGeneral/nenbar/projects/CITE/raw_files/3921_3963_4066/outs/possorted_genome_bam.bam --vcf $outFolder"/"$projectname".b38.vcf" --field GT --group-list /share/ScratchGeneral/nenbar/projects/CITE/raw_files/3921_3963_4066/outs/filtered_gene_bc_matrices/GRCh38/barcodes.tsv --out $outFolder"/"$projectname






