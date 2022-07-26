library(ShortRead)
library(GenomicRanges)


inFile="undetermined_short_S0_L001_R2_001.fastq"

reads<-readFastq(inFile)

#1. find what are the most common barcodes
bc<-as.character(reads@id)
bc<-gsub(".*:","",bc)
bcs<-sort(table(bc))

#ATCACGAT 
#22313
#GGGGGGGG
#20499
#AGGTCTAG
#18650
#ATTACTCG
#12339

gcFunction=function (x)
{
alf <- alphabetFrequency(x, as.prob = TRUE)
rowSums(alf[, c("G", "C")])
}
gc <- gcFunction(sread(reads))

gcL<-split(gc,bc)
gcLShort<-gcL[names(gcL) %in% c("ATCACGAT","AGGTCTAG","ATTACTCG","GGGGGGGG")]
sapply(gcLShort,mean)

#> sapply(gcLShort,mean)
# AGGTCTAG  ATCACGAT  ATTACTCG  GGGGGGGG 
#0.4697669 0.4718586 0.4501779 0.527988

#> sapply(gcLShort,function(x){sum(1*(x>0.8))})
#AGGTCTAG ATCACGAT ATTACTCG GGGGGGGG 
#      82       79       64     1246