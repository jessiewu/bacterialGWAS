## Code modified from the second half of Kmer_script.R (Earle, 2014) (CHW).

#######################################################################################
## prefix: prefix of the output table path
## pheno: a binary vector representing resistance phenotype
## kmerPaths: a vector of paths of the kmer files
#######################################################################################
createPhenoTable<-function(
	prefix = NULL, 
	pheno = NULL, 
	kmerPaths = NULL){
	
	out = cbind(kmerPaths, pheno)
	kmerPhenoPath<-paste0(prefix,".kmer-pheno.txt")
	message(kmerPhenoPath)
	write.table(out, file = kmerPhenoPath,
		row=F, col=F, quote=F, sep="\t")

	return(kmerPhenoPath)
	
}


#######################################################################################
## prefix: prefix of the output table path
## kmerPhenoPath: output path of the table containing kmer paths and phenotype data
#######################################################################################
gwasKmer<-function(prefix = NULL, kmerPhenoPath = NULL, gwasKmerPatternPath = NULL, minCov = NULL){
	gwasKmerOutputPath<-paste(prefix,".gwaskmer-out", sep="")
  
  kmerPheno.df = read.table(file=kmerPhenoPath, header=T, as.is=T)
	gwasKmerPatternInputFile = paste(prefix,".gwasKmerPatternInput.txt", sep="")
  write.table(kmerPheno.df[,c('filePath','phenotype')], file=gwasKmerPatternInputFile, row.names=F, col.names=F, sep="\t", quote=F)
	gwasCommand<-paste(
	  gwasKmerPatternPath,
	  gwasKmerPatternInputFile, 
		gwasKmerOutputPath,
		minCov,
		sep =" ",
		collapse=""
	)
  
  message(gwasCommand)

	## Run gwas_kmer
	system(gwasCommand)
	return(gwasKmerOutputPath)
	
}


#######################################################################################
## prefix: Prefix of the path of the plots
## results: output path of the table containing kmer paths and phenotype data
## sg: negative log10 p-values
## od: indices of the elements in sg in decresing order of the values of sg 
#######################################################################################
kmerQQPlot<-function(
	prefix = NULL,
	results = NULL,
	sg = NULL,
	od = NULL){

	results.count = length(results)
	png(paste0(prefix,"kmer_QQ.png"))
	plot.title = paste(prefix, "resistance kmer GWAS QQ plot")
	plot(x = -log10((1:results.count)/results.count), y = sg[od], 
		log = "", type = "l", 
		xlab="Null distribution of -log10(p) values",
		ylab = "Empirical distribution of -log10(p) values",
		main = plot.title
	)
	abline(0, 1, col = 2)
	dev.off()
}

#######################################################################################
## prefix: Prefix of the path of the plots
## results: output path of the table containing kmer paths and phenotype data
## sg: negative log10 p-values
## od: indices of the elements in sg in decresing order of the values of sg 
#######################################################################################
kmerEcdfPlot<-function(
	prefix = NULL,
	results = NULL,
	sg = NULL,
	od2 = NULL){
	results.count<-length(results)
	#intWidth=floor(length(results)/10)
	#nticks=seq(from=0,by=intWidth,to=nticks*500000)

	png(paste0(prefix, "kmer_ecdf.png"),width=700,height=700)
	plot(x=c(1:results.count), y = sg[od2],
		#xaxt = "n", 
    log="", type="l", 
		xlab = "Ordered kmers", ylab = "Empirical distribution of -log10(p) values",
		main = paste(prefix, "resistance kmer GWAS ECDF plot")
	)
	#axis(side=1,at=nticks)
	dev.off()
}

#######################################################################################
## prefix: Prefix of the path of the plots 
#######################################################################################
processResults<-function(prefix = NULL){
	chisqStatPath<-paste0(prefix, ".gwaskmer-out.chisqStat.txt")
	results = as.numeric(scan(chisqStatPath, what=character(0)))
	od = order(results, decreasing=TRUE)
	od2 = order(results, decreasing=FALSE)
	sg = -log10(exp(1))*pchisq(results, 1, lower=FALSE, log=TRUE)

	kmerQQPlot(prefix = prefix, results = results, sg = sg, od = od)
	kmerEcdfPlot(prefix = prefix, results = results, sg = sg, od2 = od2)


}

#######################################################################################
## kmerPaths: A list of paths to isolate specific kmer files
## phenotype: A list of phenotype corresponding to the isolates in kmerPaths
## perfix: The prefix for output files
#######################################################################################
kmerAnalysis<-function(
	kmerPaths = NULL,
	phenotype = NULL,
	prefix = NULL,
  gwasKmerPatternPath = NULL){


	## Create a table contain phenotype data and paths to the kmer files
	kmerPhenoTabPath<-createPhenoTable(
		prefix = prefix,
		pheno = phenotype,
		kmerPaths = kmerPaths
	)

	## Performs the kmer analysis
	kmerPhenoAnalysis(
		kmerPhenoTabPath = kmerPhenoTabPath,
		prefix = prefix,
		gwasKmerPatternPath = gwasKmerPatternPath
	)

}

#######################################################################################
## kmerPhenoTabPath: A list of paths to isolate specific kmer files and a list of phenotype 
## corresponding to the isolates in kmerPaths
## perfix: The prefix for output files
#######################################################################################
kmerPhenoAnalysis<-function(
	kmerPhenoTabPath = NULL,
	prefix = NULL,
	gwasKmerPatternPath = NULL,
	minCov = NULL){


	## Performs the kmer analysis
	gwasKmerOutputPath<-gwasKmer(
		prefix = prefix, 
		kmerPhenoPath = kmerPhenoTabPath,
		gwasKmerPatternPath = gwasKmerPatternPath,
		minCov = minCov
	)

	## Process results and draw plots
	processResults(prefix = prefix)

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

#######################################################################################
#######################################################################################
#######################################################################################

##################################################################################
## Help message
##################################################################################
help = paste(
  "kmerAnalysis.R Performs association tests between kmers and the phenotype.",
  "Usage: Rscript kmerAnalysis.R -dataFile kmer_pheno.txt -prefix prefix -minCov 5 -externalSofware softwarePath.txt",
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
  message(paste(args, collapse=" "))
  stop("\nIncorrect usage\n")
}

##### Get inputs from the command line #####
inputs = getCommandLineInputMatrix(args = args)

##### Extract the command line inputs #####
dataFilePath = extractInputArgument(argName = "dataFile", commandLineInputs = inputs)
prefix = extractInputArgument(argName = "prefix", commandLineInputs = inputs)
removeKmerTxt = extractInputArgument(argName = "removeKmerTxt", commandLineInputs = inputs, default = FALSE)
minCov = extractInputArgument(argName = "minCov", commandLineInputs = inputs)
externalSoftwarePaths = extractInputArgument(argName="externalSoftware", commandLineInputs = inputs)

## Get the path of external software used
externalSoftwarePaths.df = read.table(file=externalSoftwarePaths, header=T, as.is = T)
gwasKmerPatternPath = getSoftwarePath("gwasKmerPattern", externalSoftwarePaths.df)

kmerPhenoAnalysis(kmerPhenoTabPath = dataFilePath, prefix = prefix, gwasKmerPatternPath = gwasKmerPatternPath, minCov = minCov)

message(removeKmerTxt)
if(removeKmerTxt){
	kmerTxtFilePath<-paste0(prefix,".gwaskmer-out.kmer.txt")
	system(paste("rm",kmerTxtFilePath))
}
