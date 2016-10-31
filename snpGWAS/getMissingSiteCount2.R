read_reference_gz = function(ref_file) {
	unzipped<-gzfile(ref_file)
	r = scan(unzipped, what=character(0), sep="\n")
	rcat = paste(r[2:length(r)], collapse="")
	lst<-unlist(strsplit(rcat, ""))
	close(unzipped)
	return(lst)
}


get.sequence.info<-function(filepath = NULL, missing = NULL){
	#filepath = filepaths[1]
	#missing = missing.sym
	
	seq = toupper(read_reference_gz(filepath))	
	missing = toupper(missing) 
	seq.call.index<-which(!(seq %in% missing))
	seq.missing.site.counts<-length(which(seq %in% missing))	
	
	return(
		list(
			seq.call.index=seq.call.index, 
			seq.missing.site.counts=seq.missing.site.counts
		)
	)	
}


get.site.count<-function(filepath = NULL, missing = NULL){
	seq = toupper(read_reference_gz(filepath))	
	return(length(seq))
	
}


args = commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
	stop("\nIncorrect usage\n")
}

dataFile = args[1]# "/Users/jessiewu/Documents/GWAS/ecoli/treeBuildTest/fastaPaths.txt"
prefix = args[2]
totalNoCallSiteCountPath = paste(prefix, "totalNoCallSiteCount.txt", sep=".")
noCallSitesOutputPath = paste(prefix, "noCallSites.txt", sep=".")
missingSiteOuputPath = paste(prefix, "missingSites", sep=".")
#noCallSitesOutputPath = args[2] #"/Users/jessiewu/Documents/GWAS/ecoli/treeBuildTest/noCallSites.txt"
#missingSiteOuputPath = args[3]

#filepaths = "/Users/jessiewu/Documents/gwas/mtub/snpAnalysis/testMtubMissingValues.txt"
#noCallSitesOutputPath = "/Users/jessiewu/Documents/gwas/mtub/snpAnalysis/testMtubMissingNoCalls.txt"
#missingSiteOuputPath = "/Users/jessiewu/Documents/gwas/mtub/snpAnalysis/testMtubMissingCounts.txt"


data.df= read.table(file= dataFile, header=T, as.is=T)
filepaths=data.df$filePath

call.counts<-rep(0,get.site.count(filepaths[1]))

message(paste(c("Total site count:", length(call.counts)),collpase=" "))
missing.sym<-c('N','-')

missing.site.counts<-vector(length=length(filepaths))

for(i in 1:length(filepaths)){
	seq.info<-get.sequence.info(filepaths[i],missing.sym)
	call.counts[seq.info$seq.call.index] = call.counts[seq.info$seq.call.index] + 1
	
	missing.site.counts[i] = seq.info$seq.missing.site.counts
	message(paste(c("Seq",i,"has",missing.site.counts[i],"sites missing."), collapse=" "))
	
}

no.call.sites<-which(call.counts == 0)

message(paste(c("Total number of no call sites: ",length(no.call.sites)),collapse=" "))

write(paste("Total number of no call sites:", length(no.call.sites), sep=" "),
 file = totalNoCallSiteCountPath)

write.table(
	t(no.call.sites), 
	file = noCallSitesOutputPath,
	quote = FALSE,
	sep = " ",
	row.names = FALSE, 
	col.names = FALSE
)


write.table(
	missing.site.counts, 
	file = missingSiteOuputPath,
	quote = FALSE,
	sep = ",",
	row.names = FALSE, 
	col.names = FALSE
)



