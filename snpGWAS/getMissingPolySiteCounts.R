args = commandArgs(trailingOnly = TRUE)


filepathsFile = args[1]
bipinfoPath = args[2]
snpinfoPath = args[3]
output = args[4]


filepaths = scan(file=filepathsFile, what=character(0), sep="\n")
bipinfo.df = read.table(file=bipinfoPath, header=T, as.is=T)
snpinfo.df = read.table(file=snpinfoPath, header=T, as.is=T)

polySites = union(bipinfo.df$Position, snpinfo.df$Position)
missingPolySiteCounts = vector(length=length(filepaths))
for(i in 1:length(filepaths)){
	filepath = filepaths[i]
	missingSites = scan(file=filepath, what=character(0), sep="\n")
	missingPolySites = match(missingSites,polySites)
	missingPolySites = missingPolySites[!is.na(missingPolySites)]
    missingPolySiteCounts[i] = length(missingPolySites)
    if(i%%100 ==0){
    	message(i, "complete.")
    }
	
}

write(missingPolySiteCounts, file=output, ncolumns=1, sep="\n")