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
if(length(args)!=2) {
	stop("\nIncorrect usage\n")
}

filepaths = args[1]# "/Users/jessiewu/Documents/GWAS/ecoli/treeBuildTest/fastaPaths.txt"
suffix = args[2] #"/Users/jessiewu/Documents/GWAS/ecoli/treeBuildTest/noCallSites.txt"


filepaths=scan(filepaths,what=character(0),sep="\n")

call.counts<-rep(0,get.site.count(filepaths[1]))

message(paste(c("Total site count:", length(call.counts)),collpase=" "))
missing.sym<-c('N','-')

missing.site.counts<-vector(length=length(filepaths))

for(i in 1:length(filepaths)){
	filepath = filepaths[i]
	seq = toupper(read_reference_gz(filepath))	
	nSites = which(seq=='N')
	gapSites = which(seq=='-')	
	noCallSites = union(nSites,gapSites )
	outputFile = gsub(pattern = suffix, replacement="missingSites.txt", x = filepath)
	write(noCallSites, file=outputFile, ncolumns=1, sep=",")
}