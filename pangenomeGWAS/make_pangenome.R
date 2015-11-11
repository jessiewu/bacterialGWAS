# File make_pangenome.R
# Authors: Charlesworth J. and Wilson, D. J. 
#
# Copyright (C) 2015 University of Oxford
#
# This file is part of GWAS pipeline.
# See the NOTICE file distributed with this work for additional
# information regarding copyright ownership and licensing.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
#  This is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this software package; if not, write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA  02110-1301  USA

#load libraries
library(seqinr)

help = paste(
  "",
  "make_pangenome.r creates a pan-genome from denovo assemblies, using Prodigal for annotation and CD-hit for clustering",
  "Authors: Charlesworth J. & Wilson D.J.",
  "Usage: Rscript make_pangenome.r -contigFile -prefix -prodigal -similarity -coverage -externalSoftwarePaths",
  "       E.g. Rscript make_pangenome.r -contigfiles Staph_velvet.txt -prodigal yes -similarity 0.7 -coverage 0.7 -prefix Staph1 -externalSoftwarePaths soft.txt",
  "",
  "Options:",
  "-contigFile",
  " File containing list of file paths to files of assembled contigs in fasta format.",
  "-prodigal",
  " Option to run Prodigal for gene annotation (yes or no) Default=yes",
  "-similarity",
  "Similarity threshold to use for clustering open reading frames with CD-hit.Takes a value between 0.4 and 1.0.",
  "-coverage",
  "Coverage threshold to use for clustering open reading frames with CD-hit.  Takes a value between 0.4 and 1.0. Defaults to 1.0 if not specified.",
  "-externalSoftwarePaths",
  "A tab delimited file containing the name and paths of the external software used in the analysis",
  "-prefix",
  " Output file name prefix",
  "",
  sep="\n"
)

##################################################################################
## Extract paths from the input file that contains the list of external software,
## and their paths.
##################################################################################
getSoftwarePath = function(
  softwareName = NULL, 
  pathTb.df = NULL){	
  path = pathTb.df$path[which(toupper(pathTb.df$name) == toupper(softwareName))]
  checkExistence(path)
  return(path)
}

##################################################################################
## Turn command line inputs into a n x 2 matrix
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
## The first column of the matrix contains the names of the arguments, and the second
## column contains the arguments values.
##################################################################################
extractInputArgument = function(
  argName = NULL, 
  commandLineInputs = NULL,
  default = NULL){
  argIndex = which(commandLineInputs[,1] == argName)
  argIndexCount = length(argIndex)
  if(argIndexCount == 0){
    return(default)
  }else if(argIndexCount > 1){
    stop(paste(c("The argument",argName,"has been specified multiple times!"), collapse=" "))
  }else{
    return(commandLineInputs[argIndex,2])
  }
}

##################################################################################
## Format the directory path
##################################################################################
formatDir = function(dir = NULL){
  
  dirSplitted = strsplit(dir,"/")
  if(dirSplitted[length(dirSplitted)] !="/"){
    dir = paste(dir,"/",sep="")
  }
  return(dir)
  
}

##################################################################################
##Check the existence of a file.
##################################################################################
checkExistence = function(filePaths = NULL, msg = NULL){
  message(paste(filePaths))
  doesNotExist = which(!file.exists(filePaths))
  if(length(doesNotExist) > 0){
    if(!is.null(msg)){
      message(msg)
    }
    stop(paste("The file", filePaths[doesNotExist]," does not exist!", collapse="\n", sep=""))
  }
}

##################################################################################
## Get the path of a file given its directory and name.
##################################################################################
formatPath = function(dir = NULL, filename = NULL){
  dir = formatDir(dir)
  path = paste(dir,filename, sep="")
  checkExistence(path)
  return(path)
}

##################################################################################
## Save input and later run information.
##################################################################################
displayInputInfo = function(
  contigFile = NULL,
  prodigal= NULL,
  prefix= NULL,
  similarity = NULL, 
  coverage = NULL, 
  run_dir = NULL,
  inputLogFilePath = NULL){
    cat(paste0("Contig file: ", contigFile),
        paste0("Output file name prefix: ",prefix),
        paste0("Similarity: ",similarity),
        paste0("Coverage ",coverage),
        paste0("Ran in directory: ",run_dir),
        paste0("Start time: ", Sys.time()),
        file = inputLogFilePath, sep="\n"
    )
}


##################################################################################
## Rename annotated contigs so they match genome ids
##################################################################################
renameContigs=function(file=NULL){
    assembly = read.fasta(paste0(file,"_prodigal.faa"), seqtype="AA")
    contigNames = attr(assembly,"name")
    shortContigIDs= strsplit(contigNames,"_")
    genes=vector()
    for (j in 1:length(contigNames)){
      gene = paste(file,shortContigIDs[[j]][length(shortContigIDs[[j]])-1],shortContigIDs[[j]][length(shortContigIDs[[j]])],sep="_")
      genes = c(genes,gene)
    }
    write.fasta(assembly,genes,paste0(file,"_prodigal.faa"))
}


##################################################################################
## Run Prodigal https://github.com/hyattpd/Prodigal
##################################################################################
runProdigal=function(
    prodigalPath=NULL,
    data.df=NULL,
    prefix=NULL
)
  {
    for (i in 1:nrow(data.df)){
    outFile=paste0(data.df[i,]$name,"_prodigal.faa")
    system(paste0(prodigalPath," -i ",data.df[i,]$filePath," -a ",outFile," -c -m -o prodigal.out -q"))
    renameContigs(data.df[i,]$name)
    system("rm prodigal.out")
    }
    system(paste0("cat *prodigal.faa > ",prefix,"_allprot.faa"))
    system("rm *prodigal.faa")
    cat(paste0("Prodigal annotation completed: ",Sys.time()), 
        file = inputLogFilePath, append=TRUE, sep="\n")    
}


##################################################################################
##Get word size for CD-HIT
##################################################################################
getWordSize=function(similarity=NULL){
  if (similarity <= 1 && similarity > 0.7){
    wsize = 5
  }
  else if (similarity <= 0.7 && similarity > 0.6){
    wsize = 4
  }
  else if (similarity <= 0.6 && similarity > 0.5){
    wsize = 3
  }
  else if (similarity <=0.5){
    wsize= 2
  }
  return(wsize)
}  

##################################################################################
## Run CD-HIT https://github.com/weizhongli/cdhit
##################################################################################
runCdhit = function(
  cdhitPath=NULL,
  prefix=NULL,
  similarity=NULL,
  coverage=NULL,
  wsize=NULL
  )
{
  input=paste0(prefix,"_allprot.faa")
  wsize=getWordSize(similarity)
  if(!file.exists(input)){
    stop(paste0(input," file not created properly."))
  }
  system(paste0(cdhitPath," -i ",input," -o ",prefix," -c ",similarity," -d 100 -M 0 -G 1 -aL ",coverage," -aS ",coverage," -n ",wsize," -g 1"))
  system(paste0("mv ",prefix," ",prefix,".pangenome.fasta"))
  cat(paste0("Clustering with CD-HIT completed: ",Sys.time()), 
      file = inputLogFilePath, append=TRUE, sep="\n")
}



#read in command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args!=0)){
  if(args[1]=="-help" | args[1]=="-h"){
    cat(help,sep="\n")
    q("no")
  }
}

if((length(args)%%2)!=0 | length(args)==0) {
  cat(help,sep="\n")
  stop("\nNo arguments supplied\n")
}


#message("flag")
inputs = getCommandLineInputMatrix(args = args)
contigFile = extractInputArgument(argName="contigFile", commandLineInputs = inputs, default = c(0))
similarity = extractInputArgument(argName="similarity", commandLineInputs = inputs, default=c(0))
coverage= extractInputArgument(argName="coverage", commandLineInputs = inputs, default = 1.0)
prodigal= extractInputArgument(argName="prodigal", commandLineInputs = inputs, default = "yes")
prefix = extractInputArgument(argName="prefix", commandLineInputs = inputs, default = c(0))
externalSoftwarePaths = extractInputArgument(argName="externalSoftwarePaths", commandLineInputs = inputs)
run_dir= getwd()

## Read in the paths of the contig files and check that the paths are valid##
checkExistence(contigFile)
data.df = read.table(file =contigFile, header = T, as.is = T)
fastaFiles = data.df$filePath

##check option specified for prodigal input##
prodigal=tolower(prodigal)
if(prodigal!= "yes" & prodigal !="no"){
  stop("\nIncorrect format: -prodigal [yes] [no]\n")
}

##check CD-HIT parameter input values##
similarity=as.numeric(similarity)
if(similarity < 0.4 || similarity >  1){
  stop("\nIncorrect format: -similarity not in range 0.4 to 1.0\n")
}  
coverage=as.numeric(coverage)
if(coverage < 0.4 || coverage >  1){
  stop("\nIncorrect format: -coverage not in range 0.4 to 1.0\n")
}  

##### Get the path of external software used #####
externalSoftwarePaths.df = read.table(file=externalSoftwarePaths, header=T, as.is = T) 
prodigalPath = getSoftwarePath("Prodigal", externalSoftwarePaths.df)
cdhitPath = getSoftwarePath("CD-HIT", externalSoftwarePaths.df)

##### Save input and later run information #####
inputLogFilePath = paste(prefix, "input_logfile.txt", sep="_")
displayInputInfo(
  contigFile = contigFile, 
  prodigal = prodigal, 
  prefix = prefix, 
  similarity = similarity,
  coverage = coverage,
  run_dir = run_dir,
  inputLogFilePath = inputLogFilePath
)
message(paste(c("Saved input and run information to ",inputLogFilePath, "."), collapse=""))

##run Prodigal if this option selected##
if(prodigal=="yes"){
  runProdigal(prodigalPath,data.df,prefix)
}

##run CD-HIT to cluster open reading frames##
runCdhit(cdhitPath, prefix, similarity, coverage)

#read and process clusters##
clusterfile = paste0(prefix,".clstr")
checkExistence(filePaths=clusterfile, msg="The CD-HIT cluster file was not properly created")
cn <- paste("V",1:5,sep="")
clusters = read.table(clusterfile, fill=T, stringsAsFactors=F,col.names=cn)
is.newclust=clusters[,1]==">Cluster"

#get names of cluster representatives to use as ids for pan-genome
clusterrep=which(clusters[,4]=="*")
clusternames=substr(clusters[clusterrep,]$V3,1,nchar(clusters[clusterrep,]$V3)-3)
clusterid = cumsum(is.newclust)
clusterid = clusterid[clusters[,1]!=">Cluster"]
clusterid=as.factor(clusterid)
genomeid=substr(clusters[,3],2,10)
genomeid=gsub("_$","",genomeid)
genomeid=genomeid[clusters[,1]!=">Cluster"]
tb = table(genomeid,clusterid)
genomeid=as.factor(genomeid)


#prepare matrix for input to GWAS code; ##
#this is a matrix of 1s and 0s with an id column and only includes variable genes. ##
#Rows are genes in pan-genome and columns are isolate ids.##
pa = table(clusterid,genomeid)
row.names(pa)=clusternames
pa[pa>1]=1 #recode any clusters with multiple hits in a given genome to just count presence once
varGenes = pa[which(rowSums(pa)!=ncol(pa)),]
cat(ncol(varGenes))
names(varGenes)=c("geneIDs",levels(genomeid))
write.table(varGenes,paste0(prefix,".pangenome.varGenes"), row.names=T,col.names=T,quote=F,sep="\t")
cat(paste0("Pan-genome creation completed: ",Sys.time()), file=inputLogFilePath, append=TRUE, sep="\n")

