## Script to parse ncbi gene files
message("flag1")
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

help=paste("","parse_ncbi_gene.R",
           "Sarah Earle (2014)",
           "",
           "parse_ncbi_gene.R creates a BLAST database of a group of fasta files, and reads in an NCBI gene summary file to create a readable database to use to annotate BLAST results",
           "",
           "Usage: Rscript parse_ncbi_gene.R -refseq_fasta1 refseq_paths.txt -ncbi_summary genefile -prefix output_prefix","",
           "-refseq_fasta1"," List of paths to fasta files to create BLAST database with",
           "-refseq_fasta2"," List of paths to fasta files to create a second BLAST database with",
           "-ncbi_summary"," Summary file downloaded from NCBI gene (http://www.ncbi.nlm.nih.gov/gene/)",
           "-prefix"," Output file prefix",sep="\n")

args = commandArgs(trailingOnly = TRUE)

##### Get inputs from the command line #####
inputs = getCommandLineInputMatrix(args = args)

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
refseq_fasta1 = extractInputArgument(argName = "refseq_fasta1", commandLineInputs = inputs, canBeNULL = TRUE)
message(paste0("refseq_fasta1:", refseq_fasta1))
refseq_fasta2 = extractInputArgument(argName = "refseq_fasta2", commandLineInputs = inputs, canBeNULL = TRUE)
ncbi_summary = extractInputArgument(argName = "ncbi_summary", commandLineInputs = inputs)
prefix = extractInputArgument(argName = "prefix", commandLineInputs = inputs)


## Get the path of external software used
externalSoftwarePaths = extractInputArgument(argName="externalSoftware", commandLineInputs = inputs, canBeNULL = TRUE)
makeblastdbPath = NULL
if(!(is.null(refseq_fasta1) & is.null(refseq_fasta2))){
  checkExistence(externalSoftwarePaths)
  externalSoftwarePaths.df = read.table(file=externalSoftwarePaths, header=T, as.is = T)
  makeblastdbPath = getSoftwarePath("makeblastdb", externalSoftwarePaths.df)
}


checkExistence(ncbi_summary)
if(!is.null(refseq_fasta1)){
  checkExistence(refseq_fasta1)
}
if(!is.null(refseq_fasta2)){
  checkExistence(refseq_fasta2)
}


nlines=system(paste0("wc -l ",ncbi_summary),intern=TRUE);nlines=as.numeric(unlist(strsplit(nlines," "))[1])
if(nlines>1000000){
	split=floor(nlines/1000000)
} else {
	split=1
}

split=split+1

gene=c()

for(i in 1:split){
	gene=c(gene,scan(paste0(ncbi_summary),what=character(0),skip=(i-1)*1000000,nlines=1000000,sep="\n"))
}

if(length(refseq_fasta1)!=0){
existence=c()
for(i in 1:length(refseq_fasta1)){
	existence[i]=file.exists(refseq_fasta1[i])
}
if(length(which(existence==FALSE))!=0){
	w=paste(which(existence==FALSE),collapse=",")
	stop("\nIncorrect usage: Refseq fasta1 files ",w," don't exist")
}

	refseq_fasta1=scan(refseq_fasta1,what=character(0),sep="\n")
}

if(length(refseq_fasta2)!=0){

	existence=c()
	for(i in 1:length(refseq_fasta2)){
		existence[i]=file.exists(refseq_fasta2[i])
	}
	if(length(which(existence==FALSE))!=0){
		w=paste(which(existence==FALSE),collapse=",")
		stop("\nIncorrect usage: Refseq fasta2 files ",w," don't exist")
	}
	refseq_fasta2=scan(refseq_fasta2,what=character(0),sep="\n")
}

# From the 'summary' gene download, find the lines which begin with "ID:" to define separate entries

id=c()
	for(i in 1:length(gene)){
		t=unlist(strsplit(gene[i],""))
		if(length(t)>1){
			if(t[1]=="I"&t[2]=="D"&t[3]==":"){
				id[(length(id)+1)]=i
			}
		}
		#if((i%%1000)==0) cat("Done",i,"\n")
	}
	
# Output file	
genefile=paste0(prefix,".ncbigene_annotation.txt")

headerinfo=c("Index","Gene","Product","Other_aliases","Genomic_context","Accession","Pos1","Pos2","Other_designations")

cat(paste(c(headerinfo),collapse="\t"),file=genefile,append=FALSE,sep="\n")

# For each entry, parse into output file
for(i in 1:(length(id)-1)){
	
	annotation.a=c()
	annotation.b=c()
	annotation.c=c()
	ann1=c()
	ann2=c()
	ann3=c()
	ann4=c()
	ann5=c()
	ann6=c()
	
	## Pull out single record i
	rec=cbind(gene[(id[i]+1):(id[i+1])])
	# Split into space separated
	rec=unlist(strsplit(rec," "))
	if(length(which(rec=="Annotation:"))!=0){
	
		if(length(which(rec=="Aliases:"))>0){
			info=rec[3:(which(rec=="Aliases:")-2)]
			info=paste(info,collapse=" ")
			} else if(length(which(rec=="Designations:"))>0){
				info=rec[3:(which(rec=="Designations:")-2)]
				info=paste(info,collapse=" ")
				} else if(length(which(rec=="context:"))>0){
				info=rec[3:(which(rec=="context:")-2)]
				info=paste(info,collapse=" ")
				}
	
		if(length(which(rec=="Aliases:"))>0){
			if(length(which(rec=="Designations:"))>0){
				aliases=paste(rec[(which(rec=="Aliases:")+1):(which(rec=="Designations:")-2)],collapse=" ")
				} else if(length(which(rec=="context:"))>0){
					aliases=paste(rec[(which(rec=="Aliases:")+1):(which(rec=="context:")-2)],collapse=" ")
					} else if(length(which(rec=="Annotation:"))>0){
						aliases=paste(rec[(which(rec=="Aliases:")+1):(which(rec=="Annotation:")[1]-1)],collapse=" ")
						} else{
							aliases="NA"
							}
						}
	
		if(length(which(rec=="context:"))>0){
			context=paste(rec[(which(rec=="context:")+1):(which(rec=="Annotation:")[1]-1)],collapse=" ")
			} else {
				context="NA"
			}
	
		if(length(which(rec=="Annotation:"))==1){
			if(length(which(rec=="ID:"))>1){
				annotation.a=rec[(which(rec=="Annotation:")+1):(which(rec=="ID:")[length(which(rec=="ID:"))]-1)]
				} else {
					annotation.a=rec[(which(rec=="Annotation:")+1):(which(rec=="ID:")-1)]
				}
			} else if(length(which(rec=="Annotation:"))==2){
				if(length(which(rec=="ID:"))>1){
					annotation.a=rec[(which(rec=="Annotation:")[1]+1):(which(rec=="Annotation:")[2]-1)]
					annotation.b=rec[(which(rec=="Annotation:")[2]+1):(which(rec=="ID:")[length(which(rec=="ID:"))]-1)]
					} else {
						annotation.a=rec[(which(rec=="Annotation:")[1]+1):(which(rec=="Annotation:")[2]-1)]
						annotation.b=rec[(which(rec=="Annotation:")[2]+1):(which(rec=="ID:")-1)]
						}
			} else {
				if(length(which(rec=="ID:"))>1){
					annotation.a=rec[(which(rec=="Annotation:")[1]+1):(which(rec=="Annotation:")[2]-1)]
					annotation.b=rec[(which(rec=="Annotation:")[2]+1):(which(rec=="Annotation:")[3]-1)]
					annotation.c=rec[(which(rec=="Annotation:")[3]+1):(which(rec=="ID:")[length(which(rec=="ID:"))]-1)]
			
					} else {
						annotation.a=rec[(which(rec=="Annotation:")[1]+1):(which(rec=="Annotation:")[2]-1)]
						annotation.b=rec[(which(rec=="Annotation:")[2]+1):(which(rec=="Annotation:")[3]-1)]
						annotation.c=rec[(which(rec=="Annotation:")[3]+1):(which(rec=="ID:")-1)]
						}
						}

	
			if(length(which(rec=="Annotation:"))==1){
				annotation2=unlist(strsplit(annotation.a,""))
				annotation2=annotation2[(which(annotation2=="(")+1):(which(annotation2==")")-1)]
				ann.dots=which(annotation2==".")
				ann3="NA"
				ann4="NA"
				ann5="NA"
				ann6="NA"
	
			if(length(which(annotation2==","))==0){
				ann1=paste(annotation2[1:(ann.dots[1]-1)],collapse="")
				ann2=paste(annotation2[(ann.dots[2]+1):(length(annotation2))],collapse="")
				} else if(length(which(annotation2==","))!=0){
					ann1=paste(annotation2[1:(ann.dots[1]-1)],collapse="")
					ann2=paste(annotation2[(ann.dots[2]+1):(which(annotation2==",")-1)],collapse="")
					}
	
				} else if(length(which(rec=="Annotation:"))>1){
					annotation2=unlist(strsplit(annotation.a,""))
				annotation2=annotation2[(which(annotation2=="(")+1):(which(annotation2==")")-1)]
				ann.dots=which(annotation2==".")

					ann1=paste(annotation2[1:(ann.dots[1]-1)],collapse="")
				ann2=paste(annotation2[(ann.dots[2]+1):(length(annotation2))],collapse="")
				
					annotation2=unlist(strsplit(annotation.b,""))
					annotation2=annotation2[(which(annotation2=="(")+1):(which(annotation2==")")-1)]
					ann.dots=which(annotation2==".")
					ann5="NA"
					ann6="NA"
	
		if(length(which(annotation2==","))==0){
		ann3=paste(annotation2[1:(ann.dots[1]-1)],collapse="")
	ann4=paste(annotation2[(ann.dots[2]+1):(length(annotation2))],collapse="")
	} else if(length(which(annotation2==","))!=0){
		ann3=paste(annotation2[1:(ann.dots[1]-1)],collapse="")
	ann4=paste(annotation2[(ann.dots[2]+1):(which(annotation2==",")-1)],collapse="")
	
	}
	}
	
	if(length(which(rec=="Annotation:"))==3){
		annotation2=unlist(strsplit(annotation.c,""))
	annotation2=annotation2[(which(annotation2=="(")+1):(which(annotation2==")")-1)]
	ann.dots=which(annotation2==".")
	
		if(length(which(annotation2==","))==0){
		ann5=paste(annotation2[1:(ann.dots[1]-1)],collapse="")
	ann6=paste(annotation2[(ann.dots[2]+1):(length(annotation2))],collapse="")
	} else if(length(which(annotation2==","))!=0){
		ann5=paste(annotation2[1:(ann.dots[1]-1)],collapse="")
	ann6=paste(annotation2[(ann.dots[2]+1):(which(annotation2==",")-1)],collapse="")
	
	} 
	}
	
	if(length(which(rec=="Designations:"))>0){
		if(length(which(rec=="context:"))>0){
		des=paste(rec[(which(rec=="Designations:")+1):(which(rec=="context:")-2)],collapse=" ")
		} else if(length(which(rec=="context:"))==0){
			des=paste(rec[(which(rec=="Designations:")+1):(which(rec=="Annotation:")-1)],collapse=" ")
	} 
	} else {
		des="NA"
	
	}
	if(ann3=="NA"&ann4=="NA"&ann5=="NA"&ann6=="NA"){
		all.info=matrix(cbind(i,rec[2],info,aliases,context,annotation.a[1],ann1,ann2,des),ncol=9)
	} else if(ann5=="NA"&ann6=="NA"){
		all.info=rbind(all.info=matrix(cbind(i,rec[2],info,aliases,context,annotation.a[1],ann1,ann2,des),ncol=9),
			matrix(cbind(i,rec[2],info,aliases,context,annotation.a[1],ann3,ann4,des),ncol=9))
	} else {
		all.info=rbind(matrix(cbind(i,rec[2],info,aliases,context,annotation.a[1],ann1,ann2,des),ncol=9),
			matrix(cbind(i,rec[2],info,aliases,context,annotation.a[1],ann3,ann4,des),ncol=9),
			matrix(cbind(i,rec[2],info,aliases,context,annotation.a[1],ann5,ann6,des),ncol=9))
	}
	write.table(all.info,file=genefile,append=TRUE,sep="\t",quote=FALSE,row=F,col=F)
	if((i%%100)==0) cat("Done",i,"\n")
	}
}

## Fasta files from ncbi -> BLAST databases
if(length(refseq_fasta1)!=0){

	output_file=paste0(prefix,"refseq_fasta_db1_input.fa")

	for(i in 1:length(refseq_fasta1)){
		f=scan(refseq_fasta1[i],what=character(0),sep="\n")
		cat(f,file=output_file,append=TRUE,sep="\n")
	}

	## Make BLAST db1
	system(paste0(makeblastdbPath, " -in ", prefix,"refseq_fasta_db1_input.fa -dbtype nucl -title ",prefix,"_blastdb -out ",prefix,"_blastdb"))
	message("Blast databases #1 have been successfully created.")
}

if(!is.null(refseq_fasta2)){
	
	output_file=paste0(prefix,"refseq_fasta_db2_input.fa")

	for(i in 1:length(refseq_fasta2)){
		f=scan(refseq_fasta2[i],what=character(0),sep="\n")
		cat(f,file=output_file,append=TRUE,sep="\n")
	}

	## Make BLAST db2
	system(paste0(makeblastdbPath, " -in ", prefix,"refseq_fasta_db2_input.fa -dbtype nucl -title ",prefix,"_2_blastdb -out ",prefix,"_blastdb"))
	message("Blast databases #2 have been successfully created.")
}

