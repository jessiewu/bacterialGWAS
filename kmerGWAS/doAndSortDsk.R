## Code modified from do_and_sort_dsk.R (Wilson, 2014) (CHW)

# Get relative position and ordering of core SNPs in a collection of velvet assemblies

stripReads<-function(sourceFile = NULL, filePathPrefix = NULL, samtoolsPath = NULL){

	outputBam = paste0(filePathPrefix,".bam")
  system(paste(c(samtoolsPath, "view -b -F 1024 ", sourceFile, ">", outputBam), collapse=" "))
	return(outputBam)
	
}


convertBamToFastq<-function(
	filePathPrefix=NULL, 
	bamFile = NULL,
	samToFastqPath = NULL){		

	fastq1 = paste0(filePathPrefix,"_1.fastq")
	fastq2 = paste0(filePathPrefix,"_2.fastq")

	command<-paste(c(
		"java -Xmx2g -jar ", 
		samToFastqPath,
    " I=", bamFile,
		" F=", fastq1,
		" F2=", fastq2,
		" INCLUDE_NON_PF_READS=true QUIET=true VALIDATION_STRINGENCY=SILENT"), collapse="")
	system(command)

	fastqList<-list(fastq1 = fastq1, fastq2 = fastq2)

	return(fastqList)
}



trimAdaptors<-function(
	filePathPrefix = NULL,
	qualityScore = NULL,
	fastq1Path = NULL,
	fastq2Path = NULL,
  trimmomaticPath = NULL,
  trimmomaticPEPath = NULL){

	forward_paired_fastq = paste0(filePathPrefix,"_output_forward_paired.fastq")
	forward_unpaired_fastq = paste0(filePathPrefix,"_output_forward_unpaired.fastq")
	reverse_paired_fastq = paste0(filePathPrefix,"_output_reverse_paired.fastq")
	reverse_unpaired_fastq = paste0(filePathPrefix,"_output_reverse_unpaired.fastq")

	trimAdaptorCommand = NULL
	illuminaClip = paste(c("ILLUMINACLIP:", trimmomaticPEPath,":2:30:10"), collapse="")
  
	if(is.na(qualityScore)){
		trimAdaptorCommand<-paste(
			"java -jar", trimmomaticPath, trimmomaticPEPath, "-threads 1",
			fastq1Path,
			fastq2Path,
			forward_paired_fastq,
			forward_unpaired_fastq,
			reverse_paired_fastq,
			reverse_unpaired_fastq,
			illuminaClip, 
			"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
			collapse=" "
		)	
	}else{
    
    
    
		trimAdaptorCommand<-paste(
		  "java -jar", trimmomaticPath, "PE", "-threads 1",
			paste("-",qualityScore,sep=""),
			fastq1Path,
			fastq2Path,
			forward_paired_fastq,
			forward_unpaired_fastq,
			reverse_paired_fastq,
			reverse_unpaired_fastq,
			illuminaClip, 
      "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
			collapse=" "
		)
	
	}

	system(trimAdaptorCommand)

	fastqList<-list(
		forward_paired_fastq = forward_paired_fastq,
		forward_unpaired_fastq = forward_unpaired_fastq,
		reverse_paired_fastq = reverse_paired_fastq,
		reverse_unpaired_fastq = reverse_unpaired_fastq
	)

	return(fastqList)

}



interleaveFastqFiles<-function(
	filePathPrefix = NULL,
	forwardPairedFastqPath = NULL,
	forwardUnpairedFastqPath = NULL,
	reversePairedFastqPath = NULL,
	reverseUnpairedFastqPath = NULL,
  shuffleSequencesFastqPath = NULL){
	

	both_paired_fastq = paste0(filePathPrefix,"_output_both_paired.fastq")
	both_unpaired_fastq = paste0(filePathPrefix,"_output_both_unpaired.fastq")	

	pairedCommand = paste(
	  shuffleSequencesFastqPath,
		forwardPairedFastqPath,
		reversePairedFastqPath,
		both_paired_fastq,
		sep = " ",
		collapse=""
	)	

	unpairedCommand = paste(
	  shuffleSequencesFastqPath,
		forwardUnpairedFastqPath,
		reverseUnpairedFastqPath,
		both_unpaired_fastq,
		sep = " ",
		collapse=""
	)	

	system(pairedCommand)
	system(unpairedCommand)

	bothPairingFastqList<-list(
		both_paired_fastq = both_paired_fastq,
		both_unpaired_fastq = both_unpaired_fastq)

	return(bothPairingFastqList)

}



doFullDodsk<-function(filePathPrefix = NULL){
	
	full_dodsk = paste0(filePathPrefix,".solid_kmers_binary")
	full_dodsk_reads = paste0(filePathPrefix,".reads_binary")

	message(full_dodsk)
	message(full_dodsk_reads)

	kmer_files=dir("./",pattern=glob2rx("*solid_kmers_binary"))

	message(paste("#kmer files:",length(kmer_files)))

	for(i in 1:length(kmer_files)){
		message(kmer_files[i])
	}

	remaining_comid = sample(setdiff(full_dodsk, kmer_files))

	if(length(remaining_comid) != 0){
    
		cat(remaining_comid,file=paste0(filePathPrefix,"_dodsk_remaining.txt"),sep="\n")
		stop("\nThe following comids did not complete kmer counting:\n",
		paste(remaining_comid,collapse="\n"),"\nWritten incomplete comids to file")		

	}

	dodskList<-list(
		full_dodsk = full_dodsk, 
		full_dodsk_reads = full_dodsk_reads,
		remaining_comid = remaining_comid
	)

	return(dodskList)

}



createFinalKmerFiles<-function(
	kmerFileDir = NULL,
	filePathPrefix = NULL,
	full_dodsk = NULL,
  sortDskPath = NULL){

	final.kmer31.txt.gz = paste0(kmerFileDir,  filePathPrefix, ".kmer31.txt.gz")
	final.kmer31.total.txt = paste0(kmerFileDir, filePathPrefix, ".kmer31.total.txt")

	message(final.kmer31.txt.gz)
	message(final.kmer31.total.txt)
	sortDskCommand<-paste(c(sortDskPath, full_dodsk,"2>",final.kmer31.total.txt, "| gzip -c >", final.kmer31.txt.gz), collapse=" ")
  
	system.time(system(sortDskCommand))

	finalKmerList<-list(
		final.kmer31.txt.gz = final.kmer31.txt.gz,
		final.kmer31.total.txt = final.kmer31.total.txt
	)

	return(finalKmerList)
}



checkShortComid<-function(comid = NULL, short_comid = NULL){

	comidElt<-unlist(strsplit(comid,"/"))
	indexOf<-regexpr(short_comid, comidElt[length(comidElt)])[[1]]

	if(indexOf != 1){
		warning(paste("The names of the BAM file and the kmer file do not have the same prefix: ",comid, short_comid))
	}
  
}



doAndSortDsk<-function(
	comid = NULL, 
	short_comid = NULL, 
	kmerFileDir = NULL,
	qualityScore = NULL,
  samtoolsPath = NULL,
	samToFastqPath = NULL,
  trimmomaticPath = NULL,
  trimmomaticPEPath = NULL,
	shuffleSequencesFastqPath = NULL,
  sortDskPath = NULL,
  dskPath = NULL){
		

	checkShortComid(comid, short_comid)

	tmp = tempfile(pattern=paste0(short_comid,"_"),tmpdir=".")
	tmp_bam = stripReads(comid,tmp, samtoolsPath = samtoolsPath)
	tmpFastqList<-convertBamToFastq(tmp,tmp_bam, samToFastqPath = samToFastqPath)
	tmpParingFastqList<-trimAdaptors(
		filePathPrefix = tmp,
		qualityScore = qualityScore,
		fastq1Path = tmpFastqList$fastq1,
		fastq2Path = tmpFastqList$fastq2,
    trimmomaticPath = trimmomaticPath,
		trimmomaticPEPath = trimmomaticPEPath
	)
	

	tmpBothParingFastqList<-interleaveFastqFiles(
		filePathPrefix = tmp,
		forwardPairedFastqPath = tmpParingFastqList$forward_paired_fastq,
		forwardUnpairedFastqPath = tmpParingFastqList$forward_unpaired_fastq,
		reversePairedFastqPath = tmpParingFastqList$reverse_paired_fastq,
		reverseUnpairedFastqPath = tmpParingFastqList$reverse_unpaired_fastq,
		shuffleSequencesFastqPath = shuffleSequencesFastqPath
	)


	# Count kmers
	final_fastq_list = paste0(short_comid)

	system(
		paste0("printf \"",tmpBothParingFastqList$both_paired_fastq,"\\n", tmpBothParingFastqList$both_unpaired_fastq,"\\n\" > ",final_fastq_list,"; ", dskPath," ",final_fastq_list," 31")
	)
	

	# Clean up
	unlink(tmp_bam)
	lapply(tmpFastqList, unlink)
	lapply(tmpParingFastqList, unlink)
	lapply(tmpBothParingFastqList, unlink)


	# Check that counting kmers completed
	dodskList<-doFullDodsk(short_comid)

	# do_sort_dsk	
  finalKmerList<-createFinalKmerFiles(
		kmerFileDir = kmerFileDir,
		filePathPrefix = short_comid, 
		full_dodsk = dodskList$full_dodsk,
    sortDskPath = sortDskPath
	)
	
	dosortdsk_size=as.numeric(unlist(strsplit(system(paste0("ls -l ",finalKmerList$final.kmer31.txt.gz),intern=TRUE)," "))[5])>10000

	if(!dosortdsk_size){
		stop(
			"\nThe following comids did not complete kmer sorting:\n",
			paste(dodskList$remaining_comid, collapse="\n"),
			"\nWritten incomplete comids to file"
		)
	} else {
		# Clean up		
		lapply(dodskList,unlink)
	}
	
}

##################################################################################
## Turn command line inputs into a n x 2 matrix
## @args: A string vector of command line arguments.
##################################################################################
getCommandLineInputMatrix = function(args = NULL){
  argsCount = length(args)/2
  inputs = matrix(nrow=argsCount, ncol=2)
  inputs[,1] = args[c(1:argsCount)*2 -1]
  
  if(length(which(regexpr("-",inputs[,1]) != 1)) > 0){
    stop("Argument names must start with '-'! E.g. -phylogeny.")
  }
  
  inputs[,1] = gsub("-", "", inputs[,1])
  inputs[,2] = args[c(1:argsCount)*2]
  return(inputs)
  
}

##################################################################################
## Extract input arguments.
## @argName: Name of the argument to be extracted
## @commandLineInputs: A n x 2 matrix, where n is the number of command line inputs.
## @default: default value of the input argument if not provided
## The first column of the matrix contains the names of the arguments, and the second
## column contains the arguments values.
##################################################################################
extractInputArgument = function(argName = NULL, commandLineInputs = NULL,
                                default = NULL, canBeNULL = FALSE){
  
  argIndex = which(commandLineInputs[,1] == argName)
  argIndexCount = length(argIndex)
  if(argIndexCount == 0){
    if(canBeNULL | !is.null(default)){
      return(default)
    }else{
      stop(paste(c("The argument", argName, "must be specified!"), collapse=" "))
    }
  }else if(argIndexCount > 1){
    stop(paste(c("The argument", argName, "has been specified multiple times!"), collapse=" "))
  }else{
    return(commandLineInputs[argIndex,2])
  }
  
}

################################################################################################
## Format the dir appropriately.
## @dir: The directory path to be formated properly.
################################################################################################
formatDir = function(dir = NULL){
  dirLength = nchar(dir)
  if(substr(dir, dirLength, dirLength) !='/'){
    dir = paste(dir,"/",sep="")
  }
  return(dir)
}

##################################################################################
## Extract paths of the software to be used.
## @softwareName: The name of the software to be used.
## @pathTB.df: A data frame containing a list of external software and their paths
##################################################################################
getSoftwarePath = function(
  softwareName = NULL, 
  pathTb.df = NULL){  
  
  path = pathTb.df$path[which(toupper(pathTb.df$name) == toupper(softwareName))]
  
  if(length(grep("/",path)) != 0){
    checkExistence(path)
  }else{
    warning(paste(c(softwareName, "is assumed to be either in the current working directly or that its path has already been added to the PATH enviroment variable."), 
                  collapse = " "))
  }
  
  return(path)
}


################################################################################################
## Chech the existence of file paths.
## @filePath: The path of the file to be checked.
################################################################################################
checkExistence = function(filePath = NULL){  
  
  doesNotExist = which(!file.exists(filePath))
  if(length(doesNotExist) > 0){
    stop(paste(c("The following files", filePath[doesNotExist],"do not exist!"), collapse="\n"))
  }
  
}



##################################################################################
## Help message
##################################################################################
help = paste(
"do_dsk.R Count kmers for a list of bam files, first removing duplicates and trimming adaptors",
"Usage: Rscript doAndSortDsk.R -bamFilePaths bam_file_list.txt -comids comid_list.txt -externalSofware softwarePath.txt -kmerFileDir kmerFileDir -qualityScore phred33 1> log_file.txt 2> error_log_file.txt",
"or",
"Usage: Rscript doAndSortDsk.R -bamFilePaths bam_file_list.txt -comids comid_list.txt 1> log_file.txt 2> error_log_file.txt",
sep="\n")

# Read options from command line
args = commandArgs(trailingOnly = TRUE)
if(length(args!=0)){
  if(args[1]=="-help" | args[1]=="-h"){
    cat(help,sep="\n")
    q("no")
  }
}

if((length(args)%%2)!=0 | length(args)==0) {
  cat(help,sep="\n")
  stop("\nIncorrect usage\n")
}

##### Get inputs from the command line #####
inputs = getCommandLineInputMatrix(args = args)

##### Extract the command line inputs #####
bam_list = extractInputArgument(argName = "bamFilePaths", commandLineInputs = inputs)
prefix_list = extractInputArgument(argName = "ids", commandLineInputs = inputs)
externalSoftwarePaths = extractInputArgument(argName = "externalSoftware", commandLineInputs = inputs)
kmerFileDir = formatDir(extractInputArgument(argName = "kmerFileDir", commandLineInputs = inputs, default = getwd()))
qualityScore = extractInputArgument(argName = "qualityScore", commandLineInputs = inputs, default = NA)

externalSoftwarePaths.df = read.table(file=externalSoftwarePaths, header=T, as.is = T)


## Get the path of external software used
samtoolsPath = getSoftwarePath("samtools", externalSoftwarePaths.df)
samToFastqPath = getSoftwarePath("samToFastq", externalSoftwarePaths.df)
trimmomaticPath = getSoftwarePath("trimmomatic", externalSoftwarePaths.df)
pePath = getSoftwarePath("trimmomaticPE", externalSoftwarePaths.df)
shuffleSequencesFastqPath = getSoftwarePath("shuffleSequencesFastq", externalSoftwarePaths.df)
sortDskPath = getSoftwarePath("sortDsk", externalSoftwarePaths.df)
dskPath = getSoftwarePath("dsk", externalSoftwarePaths.df)

message(paste("Kmer file directory:", kmerFileDir))

# List of comids
comid = scan(bam_list,what=character(0))
short_comid = scan(prefix_list,what=character(0))

if(length(comid) !=length(short_comid)){
	message(paste("comid length: ",length(comid)))
	message(paste("prefix length: ",length(short_comid)))
	stop("Comids and prefixes must have the same number.")
}


## do_dsk
for(i in 1:length(comid)) {
	doAndSortDsk(comid[i],short_comid[i], kmerFileDir, qualityScore, samtoolsPath = samtoolsPath, 
               samToFastqPath = samToFastqPath, trimmomaticPath = trimmomaticPath, trimmomaticPEPath = pePath,
	             shuffleSequencesFastqPath = shuffleSequencesFastqPath, sortDskPath = sortDskPath, dskPath = dskPath)
}

message("Done");

