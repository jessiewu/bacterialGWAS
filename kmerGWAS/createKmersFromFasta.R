createKmersFromFasta = function(datTab = NULL, kmerFileDir = NULL, sortDskPath = NULL, dskPath = NULL){
	
	dskOuputTable = runDSK(id = datTab$id, filepaths = datTab$filePath, dskPath = dskPath)
	
	fileCount = nrow(dskOuputTable)
	for(i in 1:fileCount){
		kmerFilePathList = createFinalKmerFiles(kmerFileDir = kmerFileDir,
			filePathPrefix = dskOuputTable[i,"id"],
			full_dodsk = dskOuputTable[i,"solidKmersBinary"],
			sortDskPath = sortDskPath)
			
			
		dosortdsk_size = file.size(kmerFilePathList$final.kmer31.txt.gz)>10000

		if(!dosortdsk_size){
			stop(
				"\nThe following file with id did not complete kmer sorting:\n",
				paste(datTab$id[i], collapse="\n")
			)
		} else {
			# Clean up	
			unlink(dskOuputTable[i,c("solidKmersBinary")])
			unlink(dskOuputTable[i,c("readsBinary")])	
		}
		message("The file ", kmerFilePathList$final.kmer31.txt.gz, " has been created successfully!")
	}
	
}

runDSK = function(id = NULL, filepaths = NULL, dskPath = NULL){
	idCount = length(id)
	filepathCount = length(filepaths)
	if(idCount != filepathCount){
		stop(paste(c("The number of ids (", idCount,") ", 
		"do not match up with the number of fasta files(", filepathCount, "). "), collapse=""))
	}
	
	dskOutputPathTable = matrix(nrow = filepathCount, ncol = 3)
	colnames(dskOutputPathTable) = c("id", "readsBinary","solidKmersBinary")
	dskOutputPathTable[,"id"] = id
	for(i in 1:filepathCount){
		system(paste(c(dskPath, filepaths[i], "31"), collapse=" "))
		dskOutputPaths = getDskOuputBinaryNames(filepath = filepaths[i])
		dskOutputPathTable[i, "readsBinary"] = dskOutputPaths$readsBinary 
		dskOutputPathTable[i, "solidKmersBinary"] = dskOutputPaths$solidKmersBinary 
	}
	
	return(dskOutputPathTable)
		
}

getDskOuputBinaryNames = function(filepath = NULL){
	solidKmersBinary = NA
	readBinary = NA
	
	filenameTemp = unlist(strsplit(filepath, split="/"))
	filename = filenameTemp[length(filenameTemp)]
	filepathStrings = unlist(strsplit(filename, split=".", fixed = TRUE))
	
	if(length(filepathStrings) == 1){
		solidKmersBinary = paste(c(filepathStrings[1], ".solid_kmers_binary"), collapse="")
		readsBinary = paste(c(filepathStrings[1], ".reads_binary"), collapse="")
	}else{
		solidKmersBinary = paste(c(filepathStrings[-c(length(filepathStrings))], "solid_kmers_binary"), collapse=".")
		readsBinary = paste(c(filepathStrings[-c(length(filepathStrings))], "reads_binary"), collapse=".")
	}
	
	return(list(solidKmersBinary = solidKmersBinary, readsBinary = readsBinary))
	
	
}

createFinalKmerFiles<-function(
	kmerFileDir = NULL,
	filePathPrefix = NULL,
	full_dodsk = NULL,
  sortDskPath = NULL){

	final.kmer31.txt.gz = paste0(kmerFileDir, "/",  filePathPrefix, ".kmer31.txt.gz")
	final.kmer31.total.txt = paste0(kmerFileDir, "/", filePathPrefix, ".kmer31.total.txt")

	sortDskCommand<-paste(c(sortDskPath, full_dodsk,"2>",final.kmer31.total.txt, "| gzip -c >", final.kmer31.txt.gz), collapse=" ")
  
	system.time(system(sortDskCommand))

	finalKmerList<-list(
		final.kmer31.txt.gz = final.kmer31.txt.gz,
		final.kmer31.total.txt = final.kmer31.total.txt
	)

	return(finalKmerList)
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
"createKmersFromFasta.R Count kmers for a list of fasta files.",
"Usage: Rscript createKmersFromFasta.R -dataFile data.txt -externalSofware softwarePath.txt -kmerFileDir kmerFileDir 1> log_file.txt 2> error_log_file.txt",
"or",
"Usage: Rscript createKmersFromFasta.R -dataFile data.txt -externalSofware softwarePath.txt 1> log_file.txt 2> error_log_file.txt",
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
dataFilePath = extractInputArgument(argName = "dataFile", commandLineInputs = inputs)
externalSoftwarePaths = extractInputArgument(argName = "externalSoftware", commandLineInputs = inputs)
kmerFileDir = formatDir(extractInputArgument(argName = "kmerFileDir", commandLineInputs = inputs, default = getwd()))

externalSoftwarePaths.df = read.table(file=externalSoftwarePaths, header=T, as.is = T)

## Get the path of external software used
sortDskPath = getSoftwarePath("sortDsk", externalSoftwarePaths.df)
dskPath = getSoftwarePath("dsk", externalSoftwarePaths.df)

message("Kmer file directory: ", kmerFileDir)

dataFile.df = read.table(file=dataFilePath, header=T, as.is=T)
createKmersFromFasta(datTab = dataFile.df, kmerFileDir = kmerFileDir, sortDskPath = sortDskPath, dskPath = dskPath)

message("All kmer files have been created successfully.")

