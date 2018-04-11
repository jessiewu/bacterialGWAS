# Annotates kmers using BLAST against both a refseq database plus the nucleotide database
# Sarah Earle 2014

## kmer annotation script
help=paste("",
           "kmer_annotation.R",
           "Usage: Rscript kmer_annotation.R -chisq_results -kmer_results -present_ctrl -present_case -prefix -blastdb1 -blastdb2 -ncbi_db -signif -nproc -LMM_kmers -LMM_output",
           "       E.g. Rscript kmer_annotation.R -chisq_results gwaskmer.chisqStat.txt -kmer_results gwaskmer.kmer.txt -present_ctrl gwaskmer.nPresentCtrl.txt -present_case gwaskmer.nPresentCase.txt -prefix annotationrun1 -blastdb1 db1 -blastdb2 db2 -ncbi_db genedb.txt -nproc 8 -LMM_kmers gwasLMM_kmers_used.txt -LMM_output gwasLMM_allkmers_out.txt",
           "",
           "Options:",
           "-chisq_results",
           " Text file containing a list of the chi-squared test statistics for each kmer",
           "-kmer_results",
           " Text file containing the kmer sequences",
           "-present_ctrl",
           " File with number of control samples each kmer is present in",
           "-present_case",
           " File with number of case samples each kmer is present in",
           "-prefix",
           " File name prefix for output files",
           "-blastdb1",
           " First blast database to use to annotate kmers",
           "-blastdb2"," Second blast database to use to annotate kmers (more extensive search)",
           "-ncbi_db",
           " Gene database to annotate blast results",
           "-signif",
           " Number of most significant kmers to annotate (Default and minimum: 1000)",
           "-nproc",
           " Number of proccessors to run blast in parallel (Default = 1)",
           "-LMM_kmers",
           " Index of kmers used for LMM analysis (ends kmers_used.txt)",
           "-LMM_output",
           " Kmer LMM output for all kmers (ends allkmers_out.txt)",
           sep="\n")



#######################################################################################

## Functions

#######################################################################################

#######################################################################################
## Read in kmers
#######################################################################################
read_kmers=function(len_results,kmer_results){
	a=c(a,scan(kmer_results,
		skip=(len_results-1)*1000000,
		n=1000000,nlines=1000000,
		what=character(0)))	
}

#######################################################################################
## Split the kmers into chunks of 1000 for faster annotation
#######################################################################################
split_kmers=function(wh){
  splitseq = NULL
  splitseq2 = NULL
  if(length(wh) > 1000){
  	nsplit=floor(length(wh)/1000)
	  splitseq=seq(from=1,by=1000,to=1000*nsplit)
	  splitseq2=c(splitseq[2:length(splitseq)]-1,length(wh))
  }else{
    splitseq = 1
    splitseq2 = length(wh)
  }
	return(cbind(splitseq,splitseq2))
}	

#######################################################################################
## Construct a file containing the kmers for the BLAST search
#######################################################################################
construct_kmer_file=function(data_loop, kmers_signif, loop_start){
	info=paste0(">kmer",data_loop)
	cat(info,file=paste0(prefix,"temp_kmer.txt"),sep="\n",append=(data_loop>loop_start))
	cat(kmers_signif[data_loop],file=paste0(prefix,"temp_kmer.txt"),append=TRUE,sep="\n")
}

#######################################################################################
## Command to run BLAST
#######################################################################################	
run_BLAST=function(prefix,nproc,db, blastnPath = NULL){
	system(paste0(blastnPath, " -query ",prefix,"temp_kmer.txt -db ",db," -evalue 10 -word_size 10 -num_threads ",nproc," -outfmt '6 qseqid sseqid sacc pident length mismatch gapopen qstart qend evalue sstart send qseq' -max_target_seqs 1 -out ",prefix,"_blast.txt"))
}

#######################################################################################
## Read in the BLAST result	
#######################################################################################	
read_blast_result=function(blastfile){
	blast.search=scan(blastfile,what=character(0))
	blast.search=matrix(blast.search,ncol=13,byrow=TRUE)
}

#######################################################################################
## Check which kmers had a result in the blast search
#######################################################################################
kmer_check=function(kmer1000,blast.search){
	check=length(which(blast.search[,1]==kmer1000))
	return(check)
}

#######################################################################################
## If more than the top hit was returned, remove any lower hits
#######################################################################################
remove_multiple=function(which.multiples,blast.search,kmer1000){
	blast.multiple=which(blast.search[,1]==kmer1000[which.multiples])
	toremove=blast.multiple[2:length(blast.multiple)]
}


#######################################################################################
## Checks which kmers had a result and which rows need to be removed
## Uses functions kmer_check and remove_multiple
#######################################################################################
check_kmers_remove_duplicates=function(kmer1000,blast.search){
	kmer1000check=sapply(paste0("kmer",kmer1000),kmer_check,blast.search=blast.search,USE.NAMES=FALSE)
	which.multiples=which(kmer1000check>1)
	toremove=unlist(sapply(which.multiples,remove_multiple,blast.search=blast.search,kmer1000=paste0("kmer",kmer1000)))
	if(length(toremove)!=0){
		blast.search=blast.search[-toremove,]
	}
	return(blast.search)
}

#######################################################################################
## Get the accession numbers and positions from the BLAST results	
#######################################################################################	
get_blast_accessions=function(blast.search){
	# Accession numbers and positions
	blast.gen=blast.search[2]
	# E-value
	blast.e=as.numeric(blast.search[10])
	# Start position
	blast.r1=as.numeric(blast.search[11])
	# End position
	blast.r2=as.numeric(blast.search[12])
	# Percentage cover
	percover=as.numeric(blast.search[5])
	return(c(blast.gen,blast.e,blast.r1,blast.r2,percover))
}	

#######################################################################################
## Function to find entry with NC_ and position from BLAST
#######################################################################################

annotation_index = function(pos,ncbi.pos){
	sapply(pos, function(i) which(ncbi.pos[,1] <=i & ncbi.pos[,2]>=i))
}
annotation_index2 = function(pos,ncbi.pos){
	sapply(pos, function(i) which(ncbi.pos[,2] <=i & ncbi.pos[,1]>=i))
}

#######################################################################################
## Get the shortened accession number from the longer entry
#######################################################################################	
get_short_accession=function(blast.gen){
	blast.gen.split=unlist(strsplit(blast.gen,""))
	gen.w=which(blast.gen.split=="|")
	nc=paste(blast.gen.split[(gen.w[3]+1):(gen.w[4]-1)],collapse="")
}

#######################################################################################
## For a kmer matching to the gene database, gets the positions of the gene match from
##  the end of the accession number
## blast.gen: long accession entry from BLAST result, containing start and end
##  positions of the gene entry
#######################################################################################

get_positions=function(blast.gen){
	
	blast.gen.split=unlist(strsplit(blast.gen,""))
	gen.w=which(blast.gen.split=="|")
	nc=paste(blast.gen.split[(gen.w[3]+1):(gen.w[4]-1)],collapse="")
	
	# Does the accession entry include a c, if it does check where it is
	includes_c=length(which(blast.gen.split[(gen.w[4]+2):length(blast.gen.split)]=="c"))
	plus_2_c=blast.gen.split[(gen.w[4]+2)]!="c"
	# Are there only two positions in this entry
	two_positions=length(which(blast.gen.split[gen.w[4]:length(blast.gen.split)]=="-"))==1
	if(is.na(includes_c)==FALSE&is.na(plus_2_c)==FALSE&length(which(blast.gen.split[gen.w[4]:length(blast.gen.split)]=="-"))!=0){
	
	if(includes_c<=1 & plus_2_c==TRUE & two_positions==TRUE){
		pos=as.numeric(paste(blast.gen.split[(gen.w[4]+2):(which(blast.gen.split=="-")-1)],collapse=""))
		pos2=as.numeric(paste(blast.gen.split[(which(blast.gen.split=="-")[1]+1):length(blast.gen.split)],collapse=""))
	} else if(includes_c<=1 & plus_2_c==FALSE & two_positions==TRUE){
		pos=as.numeric(paste(blast.gen.split[(gen.w[4]+3):(which(blast.gen.split=="-")-1)],collapse=""))
		pos2=as.numeric(paste(blast.gen.split[((which(blast.gen.split=="-")+1)):length(blast.gen.split)],collapse=""))
	} else if(includes_c<=1 & plus_2_c==FALSE & two_positions==FALSE){
		pos1=as.numeric(paste(blast.gen.split[(gen.w[4]+3):(which(blast.gen.split=="-")[1]-1)],collapse=""))
		pos2=as.numeric(paste(blast.gen.split[((which(blast.gen.split=="-")[1]+1)):(gen.w[4]+(which(blast.gen.split[(gen.w[4]+2):length(blast.gen.split)]=="c")[2])-1)],collapse=""))
		pos3=as.numeric(paste(blast.gen.split[(gen.w[4]+which(blast.gen.split[(gen.w[4]+2):length(blast.gen.split)]=="c")[2]+2):(which(blast.gen.split=="-")[2]-1)],collapse=""))
		pos4=as.numeric(paste(blast.gen.split[(which(blast.gen.split=="-")[2]+1):length(blast.gen.split)],collapse=""))
				
		pos=min(pos1,pos2,pos3,pos4)
		pos2=max(pos1,pos2,pos3,pos4)
		rm(pos1)
		rm(pos3)
		rm(pos4)
		
	} else if(includes_c==0 & two_positions==FALSE){
		pos1=as.numeric(paste(blast.gen.split[(gen.w[4]+2):(which(blast.gen.split=="-")[1]-1)],collapse=""))
		pos2=as.numeric(paste(blast.gen.split[((which(blast.gen.split=="-")[1]+1)):(gen.w[4]+which(blast.gen.split[gen.w[4]:length(blast.gen.split)]==",")-2)],collapse=""))
		pos3=as.numeric(paste(blast.gen.split[(gen.w[4]+which(blast.gen.split[gen.w[4]:length(blast.gen.split)]==",")):(which(blast.gen.split=="-")[2]-1)],collapse=""))
		pos4=as.numeric(paste(blast.gen.split[(which(blast.gen.split=="-")[2]+1):length(blast.gen.split)],collapse=""))
				
		pos=min(pos1,pos2,pos3,pos4)
		pos2=max(pos1,pos2,pos3,pos4)
		rm(pos1)
		rm(pos3)
		rm(pos4)
	} else {
		pos=0
		pos2=0
	}
	
} else {
	pos=0
	pos2=0
}
return(c(pos,pos2,nc))
}


#######################################################################################
## Get the annotations from the gene database
## Uses the function annotation_index
## positions: a 4 column matrix of the start and end positions of each gene, accession
##  number and beginning position of the gene
## ncbi.gene: the gene database
#######################################################################################

get_annotations=function(positions,ncbi.gene){
	
	pos=as.numeric(positions[1])
	pos2=as.numeric(positions[2])
	nc=positions[3]
	blast.r1=as.numeric(positions[4])
	
	# If there are matches to the accession number in the gene database
	# Get the rows of the gene database which match the accession
	if(is.na(pos)==FALSE & is.na(pos2)==FALSE & length(which(ncbi.gene[,6]==nc))!=0){
		ncbi.match=ncbi.gene[which(ncbi.gene[,6]==nc),]
		# Get the positions of the genes of these matches
		if(length(ncbi.match)==9){
			ncbi.pos=matrix(as.numeric(unlist(ncbi.match[7:8])),ncol=2)
		} else {
			ncbi.pos=matrix(as.numeric(unlist(ncbi.match[,7:8])),ncol=2)
		}
			
		# Find which gene out of the selection with matching accession number based on position
		if(pos<pos2){
			annot=unlist(sapply((pos+blast.r1),annotation_index,ncbi.pos=ncbi.pos))[1]
		} else {
			annot=unlist(sapply((pos2+blast.r1),annotation_index,ncbi.pos=ncbi.pos))[1]
		}

		# If the positions match for one of the entries
		if(length(annot)!=0){
			if(annot!="integer(0)" & is.na(annot)==FALSE){
				if(length(ncbi.match)>9){
				annot2=ncbi.match[annot,]
				} else {
					annot2=ncbi.match
				}
				
			} else if(annot=="integer(0)" & is.na(annot)==FALSE){
				annot=annotation2(pos+blast.r1)[1]
				if(length(ncbi.match)>9){
				annot2=ncbi.match[annot,]
				} else {
					annot2=ncbi.match
				}
			} else if(is.na(annot)==TRUE){
				annot2=matrix("NA",9,nrow=1)
			}
		} else {
			annot2=matrix("NA",9,nrow=1)
		}
	} else {
		annot2=matrix("NA",9,nrow=1)
	}		


test_hash=unlist(strsplit(annot2[2],""))
if(length(which(test_hash=="#"))>0){
	test_hash[which(test_hash=="#")]="_"
	annot2[2]=paste(test_hash,collapse="")
}

return(annot2)
}

#######################################################################################
## If some kmers had a BLAST result and some didn't, order them by significance again
#######################################################################################	
order_BLAST_results=function(info){
	info_num=as.numeric(info[,1])
	info_sort=order(info_num)
	info=info[info_sort,]
	info=info[,2:ncol(info)]
	info=as.matrix(info,ncol=ncol(info))
}

#######################################################################################
## Function to call the other functions to annotate the kmers
#######################################################################################	
annotate_kmers=function(splitseq,prefix,nproc,blastdb1,blastdb2,wh_pos,kmers_signif,sg2,present_ctrl,present_case,kmer_output_file,blastnPath = NULL){
	data_loop=seq(from=splitseq[1],by=1,to=splitseq[2])
	
	### Annotate using the genus BLAST database and specified gene database
	# Construct file containing 1000 kmers
	sapply(data_loop,construct_kmer_file,kmers_signif=kmers_signif, loop_start = data_loop[1])
	## Run BLAST from the command line against refseq database and store wanted values
	run_BLAST(prefix,nproc,blastdb1, blastnPath = blastnPath)
	# Open BLAST results
	blastfile=paste0(prefix,"_blast.txt")
	blast.search=read_blast_result(blastfile)
	# kmer check - find which of the 1000 kmers have a result
	blast.search=check_kmers_remove_duplicates(data_loop,blast.search)
	# Which kmers had no BLAST result		
	check_noresult=which(sapply(paste0("kmer",data_loop),kmer_check,blast.search=blast.search,USE.NAMES=FALSE)==0)
	# Which kmers did have a BLAST result	
	whichrun=setdiff(1:1000,check_noresult)
	index=whichrun
	# Get BLAST accessions
	blast_acc=apply(blast.search,1,get_blast_accessions)
	blast_acc=matrix(blast_acc,ncol=5,byrow=TRUE)
	# Get positions
	positions=sapply(blast_acc[,1],get_positions)
	positions=matrix(positions,ncol=3,byrow=TRUE)
	# Get annotations for those with a BLAST result to the genus BLAST database
	annotation=apply(cbind(positions,blast_acc[,3]),1,get_annotations,ncbi.gene=ncbi.gene)
	annotation=matrix(annotation,ncol=9,byrow=TRUE)
	info1=cbind(whichrun,paste0("kmer",wh_pos[data_loop[1]+whichrun-1]),kmers_signif[data_loop[1]+whichrun-1],sg2[data_loop[1]+whichrun-1],present_ctrl[data_loop[1]+whichrun-1],present_case[data_loop[1]+whichrun-1],blast_acc,positions[,3],annotation)
	# When there is no match to the genus BLAST database
	if(length(check_noresult)>0){
		annotation=matrix("NA",15,nrow=length(check_noresult))
		info2=cbind(check_noresult,paste0("kmer",wh_pos[data_loop[1]+check_noresult-1]),kmers_signif[data_loop[1]+check_noresult-1],sg2[data_loop[1]+check_noresult-1],present_ctrl[data_loop[1]+check_noresult-1],present_case[data_loop[1]+check_noresult-1],annotation)
		
	}
	if(length(check_noresult==0)!=0){
		info.a=rbind(info1,info2)
	} else {
		info.a=info1
	}
	info.a=order_BLAST_results(info.a)
	
	#### Also run BLAST on second database, the whole nucleotide database
	# Run BLAST	
	run_BLAST(prefix,nproc,blastdb2, blastnPath = blastnPath)
	# Open BLAST results
	blastfile=paste0(prefix,"_blast.txt")
	blast.search=read_blast_result(blastfile)
	# kmer check - find which of the 1000 kmers have a result
	blast.search=check_kmers_remove_duplicates(data_loop,blast.search)
	# Which kmers had no BLAST result		
	check_noresult=which(sapply(paste0("kmer",data_loop),kmer_check,blast.search=blast.search,USE.NAMES=FALSE)==0)
	# Which kmers did have a BLAST result	
	whichrun=setdiff(1:1000,check_noresult)
	index=whichrun
	# Get BLAST accessions
	blast_acc=apply(blast.search,1,get_blast_accessions)
	blast_acc=matrix(blast_acc,ncol=5,byrow=TRUE)
	# Get the short accession number
	nc=sapply(blast_acc[,1],get_short_accession,USE.NAMES=FALSE)
	info3=cbind(whichrun,paste0("kmer",wh_pos[data_loop[1]+whichrun-1]),kmers_signif[data_loop[1]+whichrun-1],sg2[data_loop[1]+whichrun-1],present_ctrl[data_loop[1]+whichrun-1],present_case[data_loop[1]+whichrun-1],blast_acc,nc)	
	if(length(check_noresult)>0){
		annotation=matrix("NA",6,nrow=length(check_noresult))
		info4=cbind(check_noresult,paste0("kmer",wh_pos[data_loop[1]+check_noresult-1]),kmers_signif[data_loop[1]+check_noresult-1],sg2[data_loop[1]+check_noresult-1],present_ctrl[data_loop[1]+check_noresult-1],present_case[data_loop[1]+check_noresult-1],annotation)
	}
	
	if(length(check_noresult)!=0){
		info.b=rbind(info3,info4)
	} else {
		info.b=info3
	}
	info.b=order_BLAST_results(info.b)
	info.all=cbind(info.a,info.b)

	write.table(info.all,file=kmer_output_file,quote=F,row.names=F,col.names=F,sep="\t",append=TRUE)
	cat(paste0("Annotated ",data_loop[length(data_loop)]," kmers ",Sys.time()),file=paste0(prefix,"_",signif,"_annotlogfile.txt"),sep="\n",append=TRUE)
	cat("Done",data_loop[length(data_loop)],"\n")
	
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


checkDBExistence = function(check = NULL){
  dbpath=unlist(strsplit(check,"/"))
  db=dbpath[length(dbpath)]
  dbpath=paste(dbpath[1:(length(dbpath)-1)],collapse="/")
  existance=dir(dbpath,pattern=glob2rx(paste0(db,"*")))
  if(length(existance)<3){
    stop("\nIncorrect usage: ",check," database doesn't exist")
  }
  
}


#######################################################################################
## Read in command line arguments
#######################################################################################

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
chisq_results = extractInputArgument(argName = "chisq_results", commandLineInputs = inputs)
kmer_results = extractInputArgument(argName = "kmer_results", commandLineInputs = inputs)
present_ctrl = extractInputArgument(argName = "present_ctrl", commandLineInputs = inputs)
present_case = extractInputArgument(argName = "present_case", commandLineInputs = inputs)
prefix = extractInputArgument(argName = "prefix", commandLineInputs = inputs)
blastdb1 = extractInputArgument(argName = "blastdb1", commandLineInputs = inputs)
blastdb2 = extractInputArgument(argName = "blastdb2", commandLineInputs = inputs)
ncbi_db = extractInputArgument(argName = "ncbi_db", commandLineInputs = inputs)
signif = as.numeric(extractInputArgument(argName = "signif", commandLineInputs = inputs, default=1000))
nproc = as.numeric(extractInputArgument(argName = "nproc", commandLineInputs = inputs, default=1))
LMM_kmers = extractInputArgument(argName = "LMM_kmers", commandLineInputs = inputs)
LMM_output = extractInputArgument(argName = "LMM_output", commandLineInputs = inputs)
externalSoftwarePaths = extractInputArgument(argName="externalSoftware", commandLineInputs = inputs)


## Get the path of external software used
externalSoftwarePaths.df = read.table(file=externalSoftwarePaths, header=T, as.is = T)
blastnPath = getSoftwarePath("blastn", externalSoftwarePaths.df)
printOutTopXChisqPath = getSoftwarePath("PrintOutTopXChisq", externalSoftwarePaths.df)

#######################################################################################
## Check files exist and basic input checks
#######################################################################################

check_existance=c(chisq_results,kmer_results,present_ctrl,present_case,blastdb1,blastdb2,ncbi_db,LMM_kmers,LMM_output)


for(i in 1:length(check_existance)){
	check=check_existance[i]
	if(check!=blastdb1 & check!=blastdb2){
		if(file.exists(check)==FALSE){
			stop("\nIncorrect usage: ",check," file doesn't exist\n")
		}
	} else {
		dbpath=unlist(strsplit(check,"/"))
		db=dbpath[length(dbpath)]
		dbpath=paste(dbpath[1:(length(dbpath)-1)],collapse="/")
		existance=dir(dbpath,pattern=glob2rx(paste0(db,"*")))
		if(length(existance)<3){
			stop("\nIncorrect usage: ",check," database doesn't exist")
		}
	}
}


if(nproc>as.numeric(system("nproc",intern=TRUE))){
	stop("\nIncorrect usage: Number of processors specified is greater than the number of the machine")
}
if((nproc%%1)!=0 | nproc==0){
	stop("\nIncorrect usage: number of processors must be an integer\n")
}

if((signif%%1)!=0 | signif==0){
	stop("\nIncorrect usage: number of kmers to annotate must be an integer\n")
}


#######################################################################################
## Read in files
#######################################################################################

# Read in chi-squared statistics
results=as.numeric(scan(paste0(chisq_results),what=character(0),sep="\n"))

# Read in LMM kmer index
LMM_kmers=scan(LMM_kmers)
# Read in LMM output
LMM_output=read.table(LMM_output,header=T,sep="\t",quote="")
# Get -log10
sg=LMM_output$log10
# Get the order of significance
od=order(sg,decreasing=TRUE)
# Which are the top x kmers
odsignif=od[1:signif]
# sg2 is the -log10 of the chi squared statistic of the top x kmers in order
sg2=sg[odsignif]
# Get the kmer index
odsignif=LMM_kmers[odsignif]

# Read in number of present in cases and controls (and subset to those significant)
present_ctrl=as.numeric(scan(paste0(present_ctrl),what=character(0),sep="\n"))[odsignif]
present_case=as.numeric(scan(paste0(present_case),what=character(0),sep="\n"))[odsignif]






start.time = proc.time()

#
wh_pos = odsignif
wh_pos_output_path = paste(c(getwd(),"/",prefix,"Top",signif,"ChisqStat.txt"), collapse="", sep="")

message(paste(c("Writing to", wh_pos_output_path), collapse=" ", sep=""))	

write.table(wh_pos,
	file = wh_pos_output_path,
	quote = F, row.names = F, col.names = F
)

topXChisqPath = paste(c(getwd(),"/",prefix,"_top_",signif,"_chisqStat.txt"), collapse="", sep="")
topXKmerPath = paste(c(getwd(),"/",prefix,"_top_",signif,"_kmer.txt"), collapse="", sep="")
#topXKmerIndexPath = paste(c(getwd(),"/",prefix,"_top_",signif,"_kmerIndex.txt"), collapse="", sep="")

message(paste(c("Writing to", topXChisqPath), collapse=" ", sep=""))	
message(paste(c("Writing to", topXKmerPath), collapse=" ", sep=""))	

printOutTopXChisqPathSplit = unlist(strsplit(printOutTopXChisqPath, split = "/"))
javaFileName = printOutTopXChisqPathSplit[length(printOutTopXChisqPathSplit)]
printOutTopXChisqPathSplitDir = gsub(javaFileName, "", printOutTopXChisqPath)
javaFileName = gsub(".class", "", printOutTopXChisqPathSplit[length(printOutTopXChisqPathSplit)])

system(paste(c("java -cp",
               printOutTopXChisqPathSplitDir,
               javaFileName,
               wh_pos_output_path,
               chisq_results,
               kmer_results,
               topXChisqPath,
               topXKmerPath),
             collapse=" ", sep="")
       )


#system(paste(c("java -cp /home/wu/scripts/java/ PrintOutTopXChisqLMM", 
#			wh_pos_output_path, 
#			chisq_results, 
#			kmer_results, 
#			topXChisqPath, 
#			topXKmerPath,
#			topXKmerIndexPath
#		), 
#	collapse=" ", sep="")
#)




end.time = proc.time()
message(paste(end.time - start.time, collapse = " ",sep=""))

# wh is the top x chisquared statistics in order
wh=results[odsignif]


## Identify most significant kmers in order
topXKmers = scan(paste0(topXKmerPath), what=character(0), sep="\n")
#topXKmersIndex = scan(paste0(topXKmerIndexPath), what=numeric(0), sep="\n")
odsignifSort = sort(odsignif)

#odsignifSortPath = paste(c(getwd(),"/",prefix,"Top",signif,"_odsignifSort.txt"), collapse="", sep="")
#write.table(odsignifSort,
#	file = odsignifSortPath,
#	quote = F, row.names = F, col.names = F
#)

#message(length(topXKmersIndex))
#message(length(odsignifSort))

#message(class(topXKmersIndex))
#message(class(odsignifSort))

#message(paste(topXKmersIndex[1:100], collapse=",", sep=""))
#message(paste(odsignifSort[1:100], collapse=",", sep=""))

#message(identical(topXKmersIndex,odsignifSort))

newIndex = sapply(odsignif,function(index){which(odsignifSort==index)})
kmers_signif=topXKmers[newIndex]

# Split into chunks of 1000 to BLAST simultaneously - much faster
splitseq=split_kmers(wh)


# Read in the kmers in chunks - much faster than reading in all together
#a=c()
# Split up number of kmers into chunks of 1000000 to read in the kmers in chunks
#len_results=seq(from=1,by=1,to=(floor(length(results)/1000000)+1))
#a=sapply(len_results,read_kmers,kmer_results=kmer_results)
#a=unlist(a)

#if(length(a)!=length(results)){
#	stop("\nIncorrect usage: number of kmers doens't equal number of chi-squared results")
#}



#wh_pos=odsignif
# wh is the top x chisquared statistics in order
#wh=results[odsignif]
## Identify most significant kmers in order
#kmers_signif=a[odsignif]
# Split into chunks of 1000 to BLAST simultaneously - much faster
#splitseq=split_kmers(wh)

# Read in gene database
ncbi.gene=scan(paste0(ncbi_db),what=character(0),sep="\t",quote="")
ncbi.gene=matrix(ncbi.gene,ncol=9,byrow=TRUE)
#######################################################################################
# Output file name and header information
#######################################################################################

## Output file name
kmer_output_file=paste0(prefix,".out.",signif,".annotated.txt")
# Write the annotation start time to a log file
cat(paste0("Annotation start time ",Sys.time()),file=paste0(prefix,"_",signif,"_annotlogfile.txt"),sep="\n",append=TRUE)
# Write headers to file
#headerinfo=c("Kmer","Kmer_sequence","Pvalue","Present_ctrl","Present_cases","Accession_all","Evalue","Start_pos","End_pos","Length_match","Accession","Gene_db_index","Gene","Product","Other_aliases","Genomic_context","Accession","Gene_start","Gene_end","Gene_starta","Gene_endb","Gene_startc","Gene_endc","Other_designations","Kmer","Kmer_sequence","Pvalue","Present_ctrl","Present_cases","Accession_all","Evalue","Start_pos","End_pos","Length_match","Accession")
headerinfo=c("Kmer","Kmer_sequence","Pvalue","Present_ctrl","Present_cases","Accession_all","Evalue","Start_pos","End_pos","Length_match","Accession","Gene_db_index","Gene","Product","Other_aliases","Genomic_context","Accession","Gene_start","Gene_end","Other_designations","Kmer","Kmer_sequence","Pvalue","Present_ctrl","Present_cases","Accession_all","Evalue","Start_pos","End_pos","Length_match","Accession")
cat(paste(c(headerinfo),collapse="\t"),file=kmer_output_file,append=FALSE,sep="\n")

#######################################################################################
## Run annotation calling annotate_kmers function
#######################################################################################

apply(splitseq,1,annotate_kmers,prefix=prefix,nproc=nproc,blastdb1=blastdb1,blastdb2=blastdb2,wh_pos=wh_pos,kmers_signif=kmers_signif,sg2=sg2,present_ctrl=present_ctrl,present_case=present_case,kmer_output_file=kmer_output_file, blastnPath = blastnPath)

cat(paste0("Annotation of ",signif," top kmers ended: ",Sys.time()),file=paste0(prefix,"_",signif,"_annotlogfile.txt"),append=TRUE,sep="\n")