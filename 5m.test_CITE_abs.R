
sampleNames<-c("ETIMAS_S21", "ETIMAS_S22", "ETIMAS_T21", "ETIMAS_T22","4T1_3x")
sampleNames=c( "SVH130918A_control", "3921_3963_4066", "4497-1", "4497-2", "4290B", "4515", "4530-N", "4530-T", "MIA", "pool1", "TONS_PBMC" )

alltags<-read.table("tags_mouse.csv",sep=",")
for(sampleName in sampleNames){

	cat(sampleName)
	cat("\n")
	ABs<-read.table(paste0(sampleName,"/CITE_fastqs/tags.csv"),sep=",")
	cat(sum(ABs$V1 %in% alltags$V1)/dim(ABs)[1])
	cat("\n")
	cat(as.character(ABs[!(ABs$V1 %in% alltags$V1),"V2"]))
	cat("\n")

}