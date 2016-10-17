# kmer_LMM.R
# Sarah Earle, Chieh-Hsi Wu (2015)
# Runs an LMM analysis using a predefined relatedness matrix on kmers


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
  
  pathIndex = which(toupper(pathTb.df$name) == toupper(softwareName))
  if(length(pathIndex) == 0){
    stop(paste(softwareName, "is not specified in the external software path file.!", sep = " "))
  }
  
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
  if(is.null(filePath)){
    stop("The file path is null!")
  }
  
  doesNotExist = which(!file.exists(filePath))
  if(length(doesNotExist) > 0){
    stop(paste(c("The following files", filePath[doesNotExist],"do not exist!"), collapse="\n"))
  }
  
}

help=paste("",
           "kmer_LMM.R",
           "Usage: Rscript kmer_LMM.R chisq_results patterns pattern_index number_signif rel_matrix phenotype out_prefix",
           "",
           "chisq_results:",
           " Text file containing a list of the chi-squared statistics for each kmer",
           "patterns:",
           " Text file containing list of unique patterns of the kmers in the dataset",
           "pattern_index:",
           " Text file containing the index for each kmer of which pattern it is",
           "number_signif:",
           " Number of most significant kmers to run LMM (will be matched by an equal number of kmers randomly sampled from those remaining)",
           "rel_matrix:",
           " GEMMA relatedness matrix created using SNP data",
           "phenotype:",
           " Text file containing a list of binary phenotypes",
           "out_prefix:",
           " Output file prefix",
           sep="\n")

message("Running kmerLMM.R")
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

message(paste0("#args: ", length(args)))

##### Get inputs from the command line #####
inputs = getCommandLineInputMatrix(args = args)

##### Extract the command line inputs #####
chisq_results = extractInputArgument(argName = "chisqStat", commandLineInputs = inputs)
patterns = extractInputArgument(argName = "patternKey", commandLineInputs = inputs)
index = extractInputArgument(argName = "patternIndex", commandLineInputs = inputs)
signif = as.numeric(extractInputArgument(argName = "signif", commandLineInputs = inputs))
rel_matrix = extractInputArgument(argName = "relateMatrix", commandLineInputs = inputs)
phenotype_file = extractInputArgument(argName = "phenotype", commandLineInputs = inputs)
prefix = extractInputArgument(argName = "prefix", commandLineInputs = inputs)
externalSoftwarePaths = extractInputArgument(argName="externalSoftware", commandLineInputs = inputs)

## Get the path of external software used
externalSoftwarePaths.df = read.table(file=externalSoftwarePaths, header=T, as.is = T)
gemmaPath = getSoftwarePath("GEMMA", externalSoftwarePaths.df)



checkExistence(c(chisq_results,patterns,index,rel_matrix,phenotype_file))
checkExistence(gemmaPath)



if((signif%%1)!=0 | signif==0){
	stop("\nIncorrect usage: number of kmers to run must be an integer\n")
}



#######################################################################################
## Read in files
#######################################################################################

# Read in patterns
patterns=scan(patterns,what=character(0),sep="\n")

# Read in chi-squared statistics
results=as.numeric(scan(paste0(chisq_results),what=character(0),sep="\n"))
# Get -log10 of the chisquared statistic
sg=-log10(exp(1))*pchisq(results,1,lower=FALSE,log=TRUE)
# Get the number of samples in the dataset
pat_length=length(unlist(strsplit(patterns[1],"")))
# Read in pattern index and convert to 1 based numbers
index=scan(index,what=numeric(0),sep=",")
index=index+1
# Order kmers by significance
od=order(sg,decreasing=TRUE)
# Which are the top x kmers
odsignif=od[1:signif]


# Of the kmers not in the top x most significant, match the same number out of those remaining
rem=setdiff(od,odsignif)
message("signif:")
message(signif)

if(length(rem) > 0){
	#Get extras
	extras=sample(rem,signif,replace=FALSE)
	
	# All kmers to run through LMM
	LMM_kmers=c(odsignif,extras)
	
}else{
	LMM_kmers=odsignif
}
# Which patterns do the kmers correspond to
patterns_LMM=index[LMM_kmers]
# Reduce down to the unique patterns to run LMM on
which_patterns=unique(patterns_LMM)

#######################################################################################
## Produce GEMMA formatted files
#######################################################################################

# Produce a matrix to fill to create the GEMMA gen file
message("Produce a matrix to fill to create the GEMMA gen file:")
all_patterns=matrix(0,length(which_patterns),(pat_length+3))
for(i in 1:length(which_patterns)){
	pat=as.numeric(unlist(strsplit(patterns[which_patterns[i]],"")))
	countallele0 = length(which(pat==0))
	countallele1 = length(which(pat==1))
	
	if(countallele0<=countallele1){
		mallele=0
		nonmallele=1
	} else {
		mallele=1
		nonmallele=0
	}
	
	all_patterns[i,]=c(paste0("rs",i),mallele,nonmallele,pat)
}

# Create GEMMA snp file
message("Create GEMMA snp file:")
snp_info=cbind(paste0("rs",1:nrow(all_patterns)),c(1:nrow(all_patterns)),rep(24,length(nrow(all_patterns))))
# Find the patterns which each kmer corresponds to out of those used in the LMM
message("# Find the patterns which each kmer corresponds to out of those used in the LMM:")
LMM_pattern_index=rep(0,length(patterns_LMM))
for(i in 1:length(which_patterns)){
	LMM_pattern_index[which(patterns_LMM==which_patterns[i])]=i
}

#######################################################################################
## Output files and indexes
#######################################################################################

# Output the index to find for all kmers the pattern they correspond to
cat(LMM_pattern_index,file=paste0(prefix,"_LMM_kmer_pattern_index.txt"),sep="\n")
# Output which kmers were included in the analysis
cat(LMM_kmers,file=paste0(prefix,"_kmers_used.txt"),sep="\n")


# Output snp file
kmerGemmaSnpFormatFile = paste0(prefix,"_kmer_patterns_gemma_snp_format.txt")
write.table(snp_info,file=kmerGemmaSnpFormatFile,row=F,col=F,sep="\t",quote=FALSE)

# Output gen file
kmerGemmaGenFormatFile = paste0(prefix,"_kmer_patterns_gemma_gen_format.txt")
write.table(all_patterns,file=kmerGemmaGenFormatFile,row=F,col=F,sep="\t",quote=F)

#######################################################################################
## Run LMM
#######################################################################################


# Run LMM
message(paste0("gemmaPath: ",gemmaPath))
kmerGemmaCommand = paste0(gemmaPath, 
                          " -g ", kmerGemmaGenFormatFile, 
                          " -p ",phenotype_file,
                          " -a ", kmerGemmaSnpFormatFile, 
                          " -k ",rel_matrix," -lmm 4",
                          " -o ",prefix,"_kmer_lmmout -maf 0")
message(kmerGemmaCommand)
system(kmerGemmaCommand)


#######################################################################################
## Read in LMM output and produce summary plots
#######################################################################################
gem=read.table(paste0("./output/",prefix,"_kmer_lmmout.assoc.txt"),header=TRUE,as.is=TRUE)

gem_all=gem[LMM_pattern_index,]
log10=rep(0,nrow(gem_all))
log10[which(gem_all$p_lrt!=0)]=-log10(gem_all$p_lrt[which(gem_all$p_lrt!=0)])
gem_all=cbind(LMM_kmers,gem_all,log10)

write.table(gem_all,file=paste0(prefix,"_LMM_allkmers_out.txt"),row=F,col=T,sep="\t",quote=FALSE)

png(paste0(prefix,"_kmer_LMM_plot.png"),width=1000,height=750)
plot(x=gem_all$LMM_kmers,y=gem_all$log10,type="p",xlim=c(0,max(gem_all$LMM_kmers)),ylim=c(0,max(gem_all$log10)),xlab="Alphabetically ordered kmers",ylab="-log10(p-value)",main=paste0(prefix," Kmer LMM plot"),col=c(rep("red",length=signif),rep("blue",length=signif)))
dev.off()

png(paste0(prefix,"_kmer_chisq_lmm_comparison.png"),width=800,height=800)
plot(x=sg[LMM_kmers],y=gem_all$log10,type="p",xlab="Chi-squared -log10(p-value)",ylab="LMM -log10(p-value)",main=paste0(prefix," chi-squared against LMM"),col=c(rep("red",length=signif),rep("blue",length=signif)))
dev.off()