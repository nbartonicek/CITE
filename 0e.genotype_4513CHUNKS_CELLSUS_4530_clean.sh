#!/bin/bash

module load gi/gcc/4.8.2
module load gi/picard-tools/1.138
module load gi/plink/prebuilt/1.07
module load gi/plink/prebuilt/1.90beta_3g
module load gi/novosort/precompiled/1.03.08
module load phuluu/vcftools/0.1.15
module load gi/bedtools/2.22.0
#raw files directory
homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/CITE"

scriptsPath="$homedir/projects/CITE/scripts"
logDir=$projectDir"/scripts/logs"
mkdir -p $logDir

projectname="4513_4530_Blood"

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

(echo cel_files; grep $step1Dir/"cel_list1.txt" -v -f $step1Dir/"failed.txt" | grep "Blood" | grep '4513\|4530') >$step1Dir/"cel_list2.txt"

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
cp -R /share/ScratchGeneral/nenbar/projects/CITE/genotyping/$projectname/step1/AxiomAnalysisSuiteData /share/ScratchGeneral/nenbar/projects/CITE/genotyping/$projectname/step1/OTV
/share/ScratchGeneral/nenbar/local/lib/apt-2.10.2.2-x86_64-intel-linux/bin/apt-format-result --batch-folder $step1Dir/OTV --snp-list-file $step1Dir/OTV/OTV.keep.ps --annotation-file $annotationDir/Axiom_UKB_WCSG.na35.annot.db --export-vcf-file $outFolder"/"$projectname".b37.vcf"
#this plink is broken
#/share/ScratchGeneral/nenbar/local/lib/apt-2.10.2.2-x86_64-intel-linux/bin/apt-format-result --batch-folder $step1Dir/OTV --snp-list-file $step1Dir/OTV/OTV.keep.ps --annotation-file $annotationDir/Axiom_UKB_WCSG.na35.annot.db --export-plink-file $outFolder"/"$projectname".b37"
java -jar /share/ClusterShare/software/contrib/gi/picard-tools/1.138/picard.jar CreateSequenceDictionary R=/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg38/hg38.fa O=/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg38/hg38.dict
more $outFolder"/"$projectname".b37.vcf" | grep -v UNKNOWNPOSITION >$outFolder"/"$projectname".b37.clean.vcf"
java -jar ~/local/bin/picard.jar LiftoverVcf I=$outFolder"/"$projectname".b37.clean.vcf" O=$outFolder"/"$projectname".hg38.vcf" CHAIN=$annotationDir"/b37ToHg38.over.chain" REJECT=$annotationDir"/rejected_variants.vcf" R=/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg38/hg38.fa

awk '{gsub(/^chr/,""); print}' $outFolder"/"$projectname".hg38.vcf" > $outFolder"/"$projectname".b38.vcf"
############### with or without splink



############# clean SNPs by Seyhun ##################


INPUT=4513_4530_Blood.b37
OUTPUT=snp_missingness
mkdir $OUTPUT
SNP_RATE=0.03

# Step 1: Remove SNPs with high missingness

vcftools --vcf 3921_3963_4066_Blood.b37.clean.vcf --plink --out $INPUT
plink --file $INPUT --make-bed --out $outFolder"/"$projectname".b37" --noweb
plink --bfile $INPUT --make-bed --geno $SNP_RATE --out $OUTPUT/$OUTPUT --noweb

# Step 2 : Remove individuals with missing high number of SNPs

INPUT2=snp_missingness/snp_missingness
OUTPUT2=indiv_missingness
mkdir $OUTPUT2
MIND=0.03

plink --bfile ${INPUT2} --make-bed --mind ${MIND} --out ${OUTPUT2}/${OUTPUT2} --noweb


# Step3: Check all the individuals are labelled with corresponding sex correctly
#
#INPUT3=indiv_missingness/indiv_missingness
#OUTPUT=check_sex
#mkdir $OUTPUT3
#
#plink --bfile ${INPUT3} --check-sex --out ${OUTPUT3}/${OUTPUT3} --noweb


# Step 4: Subset chr 1 to 22 (it is hard to impute chr X, I also removed the mitochondrial ones but possible to keep them)

INPUT4=snp_missingness/snp_missingness
OUTPUT4=subset_chr1to22
mkdir $OUTPUT4
plink --bfile ${INPUT4} --chr 1-22 X --make-bed --out ${OUTPUT4}/${OUTPUT4} --allow-extra-chr --noweb


# Step 5: Remove duplicated SNPs

INPUT5=subset_chr1to22/subset_chr1to22
OUTPUT5=dedup
mkdir $OUTPUT5

plink --bfile ${INPUT5} --list-duplicate-vars --allow-extra-chr --out ${OUTPUT5}/${OUTPUT5}
plink --bfile ${INPUT5} --exclude ${OUTPUT5}/${OUTPUT5}.dupvar --make-bed --out ${OUTPUT5}/${OUTPUT5}

# Step 6: Make sure that rs ids are correct

INPUT6=dedup/dedup
OUTPUT6=rsIDconversion
RSIDFILE=XXX
mkdir $OUTPUT6

plink --bfile ${INPUT6} --update-name ${RSIDFILE} --make-bed --out ${OUTPUT6}/${OUTPUT6}


# Step 7: Chck the files against HRC reference panel (can use other panels too, check Michigan Server website for more info and download the relative files). You will use the created files in this output folder in step 8-10.

INPUT7=subset_chr1to22/subset_chr1to22
OUTPUT7=hrc_check
REFERENCE=/share/ScratchGeneral/nenbar/projects/CITE/genotyping/annotation/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
mkdir $OUTPUT7

plink --bfile ${INPUT7}  --freq --out ${OUTPUT7}/${OUTPUT7}
cp hrc_check/hrc_check.frq hrc_check/hrc_check.afreq
perl /share/ScratchGeneral/nenbar/projects/CITE/genotyping/annotation/HRC-1000G-check-bim-v4.2.pl -b ${INPUT7}.bim -f ${OUTPUT7}/${OUTPUT7}.afreq -r ${REFERENCE} -h 


# Step 8: Exclude SNPs that are in refernece panel but not in your dataset

INPUT8=subset_chr1to22/subset_chr1to22
OUTPUT8=snp2exclude
mkdir ${OUTPUT8}

plink --bfile ${INPUT8}  --exclude Exclude-subset_chr1to22-HRC.txt --make-bed --out ${OUTPUT8}/${OUTPUT8}


# Step 9: Make sure the all the snps are in the forward strand

INPUT9=snp2exclude/snp2exclude
OUTPUT9=forward_strand
REFERENCE=/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg19/hg19.fa
mkdir ${OUTPUT9}


# Get diagnostics
snpflip -b ${INPUT9}.bim -f ${REFERENCE} -o ${OUTPUT9}/${OUTPUT9}_snpflip

# Remove ambiguous SNPs
plink --bfile ${INPUT9} --exclude ${OUTPUT9}/${OUTPUT9}_snpflip.ambiguous --make-bed --out ${OUTPUT9}/${OUTPUT9}_noambiguous --noweb


# Flip SNPs with reverse strand genotypes
#plink --bfile ${OUTPUT9}/${OUTPUT9}_noambiguous --flip Strand-Flip-onek1k_rs_ids_correct-HRC.txt --make-bed --out ${OUTPUT7}/${OUTPUT7} --noweb

# Step 10: Make sure the reference alleles are correctly set

INPUT10=forward_strand/forward_strand
OUTPUT10=setREFallele
REFERENCE=/path/to/reference_file
mkdir ${OUTPUT10}

plink --bfile ${INPUT10}"_noambiguous" --a2-allele Force-Allele1-onek1k_allele_correct-HRC.txt --make-bed --out ${OUTPUT10}/${OUTPUT10}

plink --bfile ${INPUT10}"_noambiguous"  --recode vcf --out ${OUTPUT10}/${OUTPUT10}



for i in {1..22}; do 
	vcftools  --vcf  ${OUTPUT10}/${OUTPUT10}".vcf"  --chr $i  --recode --recode-INFO-all --out  ${OUTPUT10}"/"VCF_$i;
	vcf-sort ${OUTPUT10}"/VCF_"$i".recode.vcf" | /share/ScratchGeneral/nenbar/local/lib/demuxlet/bgzip -c > ${OUTPUT10}"/"$projectname".b37.$i.vcf.gz"
done





############# clean SNPs ##################

vcftools --gzvcf snps.vcf.gz --maf 0.01 --max-missing 0.8 --min-meanDP 10 --recode --out snp_filtered.vcf.gz
java -Xmx24g -jar /programs/beagle41/beagle41.jar gt=/workdir/moshood/snps/snp_filtered.vcf.gz.recode.vcf nthreads=8 niterations=0  out=imputed_file

https://www.protocols.io/view/genotype-imputation-workflow-v3-0-nmndc5e?version_warning=no

#https://github.com/statgen/demuxlet/tree/master/tutorial
vcftools --vcf 3921_3963_4066_Blood.b37.vcf --plink --out 3921_3963_4066_Blood.b37
plink --file $outFolder"/"$projectname".b37" --make-bed --out $outFolder"/"$projectname".b37" --noweb

#doesnt work liftOverPlink.py -m 3921_3963_4066_Blood.b37.map -p 3921_3963_4066_Blood.b37.ped -o 3921_3963_4066_Blood.b37.out -c $annotationDir"/b37ToHg38.over.chain"
#(1) individual and SNP missingness - per SNP/individual - SNPs missing for an individual indication of poor DNA quality, exclude SNPs missing in large proportion, usually 20%, then 0.02, first on SNPs, 
#(2) inconsistencies in assigned and genetic sex of subjects - sample mixup or procedure fail
#(3) minor allele frequency (MAF) - low allele frequence in GWAS studies a problem, but not for us
#(4) deviations from Hardy–Weinberg equilibrium (HWE) - is the genotype (AT) different from frequency of alleles in population? It should be the same unless breast cancer SNPS. HWE p value <1e−10 in cases and <1e−6 in controls.
#(5) heterozygosity rate - per individual - sample quality, 3SD
#(6) relatedness - per cohort - should be unrelated. Use pruning on autosomals
#(7) ethnic outliers (see population stratification). - per cohort


qsub -q short.q -b y -j y -N mergeBams -cwd -pe smp 20 -P TumourProgression -

type="CHUNKS"
inFile1="/paella/TumourProgressionGroupTemp/projects/data/cellranger_count/cryopreservation_paper/CID4513$type/count_CID4513_"$type"_GRCh38/outs/possorted_genome_bam.bam"
inFile2="/paella/TumourProgressionGroupTemp/projects/data/cellranger_count/human_breast/CID4530/count_CID4530_GRCh38/outs/possorted_genome_bam.bam"
outFile="/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4513"$type"_4530/outs/possorted_genome_bam.bam"
mkdir -p "/share/ScratchGeneral/nenbar/projects/CITE/raw_files/4513"$type"_4530/outs/"
samtoolsLine="samtools merge $outFile $inFile1 $inFile2"
qsub -q short.q -b y -j y -N mergeBams -cwd -pe smp 20 -P TumourProgression -V $samtoolsLine

############ submit to imputation server short

vcf-sort $outFolder"/"$projectname".b38.vcf" | /share/ScratchGeneral/nenbar/local/lib/demuxlet/bgzip -c > $outFolder"/"$projectname".b38.vcf.gz"

projectname="3921_3963_4066_Blood"
outFolder="/share/ScratchGeneral/nenbar/projects/CITE/genotyping/$projectname/step2"

for i in {1..22}; do 
	vcftools  --vcf  $outFolder"/"$projectname".hg38.vcf"  --chr chr$i  --recode --recode-INFO-all --out  $outFolder"/"VCF_$i;
	vcf-sort $outFolder"/VCF_"$i".recode.vcf" | /share/ScratchGeneral/nenbar/local/lib/demuxlet/bgzip -c > $outFolder"/"$projectname".hg38.chr$i.vcf.gz"
done

#upload to michigan imputation server, impute and download
#unzip files with password, merge, add the file used for imputation and sort
for i in {1..22}; do unzip -P P8mSMcBf5KU3ys chr_$i.zip ; done
for i in {1..22}; do unzip -P Fw4.r7eFRUgzN9 chr_$i.zip ; done
gunzip *.dose.vcf.gz
grep '^#' chr1.dose.vcf > merge.vcf
grep -v '^#' chr{1..22}.dose.vcf | sed s/^.*.vcf:// >> merge.vcf 
#grep -v '^#' $projectname.b37.vcf | sed s/^.*.vcf:// >> merge.vcf 
vcf-sort merge.vcf > merge.imputed.original.sorted.vcf

java -jar ~/local/bin/picard.jar LiftoverVcf I=$outFolder"/merge.imputed.original.sorted.vcf" O=$outFolder"/"$projectname".hg38.vcf" CHAIN=$annotationDir"/b37ToHg38.over.chain" REJECT=$annotationDir"/rejected_variants.vcf" R=/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg38/hg38.fa
java -jar ~/local/bin/picard.jar CreateSequenceDictionary R=/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg38_85/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa O=/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg38_85/Homo_sapiens.GRCh38.dna_sm.primary_assembly.dict
java -jar ~/local/bin/picard.jar LiftoverVcf I=$outFolder"/"$projectname".hg38.vcf" O=$outFolder"/"$projectname".b38.vcf" CHAIN=$annotationDir"/hg38ToGRCm38B.over.chain" REJECT=$annotationDir"/rejected_variants.vcf" R=/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg38_85/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
awk '{gsub(/^chr/,""); print}' $outFolder"/"$projectname".hg38.vcf" > $outFolder"/"$projectname".b38.vcf"

grep '^#' $outFolder"/"$projectname".b38.vcf" > $outFolder"/merge_overlapped_fullheader.vcf"

bedtools intersect -a $outFolder"/"$projectname".b38.vcf" -b /share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg38_85/Homo_sapiens.GRCh38.85.annotated.polyA.gtf >>$outFolder"/merge_overlapped_fullheader.vcf" &

#overlap with the transcriptome gtf
module load gi/bedtools/2.22.0
bedtools intersect -a $outFolder"/merge.vcf" -b /share/ClusterShare/biodata/contrib/nenbar/genomes/star/hg38_gencode.v27.polyA/gencode.v27.polyA.annotation.gtf >$outFolder"/merge_overlapped.vcf"


java -jar /share/ClusterShare/software/contrib/gi/picard-tools/1.138/picard.jar GatherVcfs I=chr{1..22}.dose.vcf O=mergedPicard.vcf


#qsub -q short.q -b y -j y -N genotypeStep4 -wd $logDir -pe smp 20 -P TumourProgression -V demuxlet --alpha 0.5 --sam /share/ScratchGeneral/nenbar/projects/CITE/raw_files/3921_3963_4066/outs/possorted_genome_bam.bam --vcf $outFolder"/"$projectname".b38.vcf" --field GT --group-list /share/ScratchGeneral/nenbar/projects/CITE/raw_files/3921_3963_4066/outs/raw_gene_bc_matrices/GRCh38/barcodes.tsv --out $outFolder"/"$projectname

#make barcodes
cat /paella/TumourProgressionGroupTemp/projects/data/cellranger_count/human_breast/CID4513/count_CID4513_GRCh38/outs/filtered_gene_bc_matrices/GRCh38/barcodes.tsv  /paella/TumourProgressionGroupTemp/projects/data/cellranger_count/human_breast/CID4530/count_CID4530_GRCh38/outs/filtered_gene_bc_matrices/GRCh38/barcodes.tsv >barcodes.tsv

sortname_line="novosort -m 16G -c 16 /share/ScratchGeneral/nenbar/projects/CITE/raw_files/genotype/test/merged.bam >/share/ScratchGeneral/nenbar/projects/CITE/raw_files/genotype/test/merged.sorted.bam"
qsub -q short.q -b y -j y -N novosort -wd $logDir -pe smp 20 -P TumourProgression -V $sortname_line

#############

############# demuxlet ##################
projectname="4513_4530_Blood"
outFolder="/share/ScratchGeneral/nenbar/projects/CITE/genotyping/$projectname/step2"

qsub -q short.q -b y -j y -N demuxlet -wd $logDir -pe smp 20 -P TumourProgression -V demuxlet --alpha 0.5 --sam /share/ScratchGeneral/nenbar/projects/CITE/raw_files/genotype/test/merged.sorted.bam --vcf $outFolder"/"$projectname".b38.vcf" --field GT --group-list /share/ScratchGeneral/nenbar/projects/CITE/raw_files/genotype/test/barcodes.tsv --out $outFolder"/"$projectname
qsub -q short.q -b y -j y -N demuxlet -wd $logDir -pe smp 20 -V demuxlet --alpha 0.5 --sam /share/ScratchGeneral/nenbar/projects/CITE/raw_files/genotype/test/merged.sorted.bam --vcf $outFolder"/merge_overlapped_fullheader.vcf" --field GT --group-list /share/ScratchGeneral/nenbar/projects/CITE/raw_files/genotype/test/barcodes.tsv --out $outFolder"/"$projectname"_overlapped"

projectname="3921_3963_4066_Blood"
outFolder="/share/ScratchGeneral/nenbar/projects/CITE/genotyping/$projectname/step2"

qsub -q short.q -b y -j y -N demuxlet -wd $logDir -pe smp 20 -V demuxlet --alpha 0.5 --sam /share/ScratchGeneral/nenbar/projects/CITE/raw_files/3921_3963_4066/outs/possorted_genome_bam.bam --vcf $outFolder"/merge_overlapped.vcf" --field GT --group-list /share/ScratchGeneral/nenbar/projects/CITE/raw_files/3921_3963_4066/outs/raw_gene_bc_matrices/GRCh38/barcodes.tsv --out $outFolder"/"$projectname"_overlapped"



