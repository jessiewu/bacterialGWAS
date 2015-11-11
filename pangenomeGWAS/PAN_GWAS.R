# File PAN_GWAS.R
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


help = paste(
  "",
  "PAN_GWAS.R runs basic and optional lmm GWAS using presence or absence of each gene in a pan-genome as the genotype",
  "Authors: Charlesworth J. & Wilson D.J.",
  "Usage: Rscript PAN_GWAS.R, -pangenome  -genomes -prefix -externalSoftwarePaths -gemma -relm",
  "       E.g. Rscript PAN_GWAS.R -contigfiles Staph_velvet.txt -similarity 0.7 -coverage 0.7 -prefix Staph1 -externalSoftwarePaths soft.txt",
  "","Options:",
  "-pangenome",
  "Pangenome matrix containing column of gene ids, column for each isolate encoding the presence or absence of each gene as 1=present or 0=absent",
  "-phenotype",
  "Matrix of phenotypes. First column contains genome ids followed by columns of binary phenotypes coded as 0=control and 1=case",
  "-gemma",
  "Run Gemma to control for population structure. (yes or no) Default = no)",
  "-relm",
  "Gemma relatedness matrix calculated from reference-based SNPs called for same isolates. Required if running Gemma",
  "-prefix",
  " Output file name prefix",
  "-externalSoftwarePaths",
  "A tab delimited file containing the name and paths of the external software used in the analysis",
  "-scriptDir",
  " Absolute path of the directory containing the scripts (Default: Current working directory)",
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
  pangenome = NULL,
  phenotype= NULL,
  prefix= NULL,
  gemma = NULL, 
  relm = NULL,
  scriptDir = NULL,
  runDir = NULL,
  inputLogFilePath = NULL){
  cat(paste0("Pangenome file: ", pangenome),
      paste0("Phenotype file: ", phenotype),
      paste0("Output file name prefix: ",prefix),
      paste0("Gemma: ",gemma), 
      paste0("Relatedness matrix file: ",relm), 
      paste0("Scripts in directory: ",scriptDir),     
      paste0("Ran in directory: ",runDir),
      paste0("Start time: ", Sys.time()),
      file = inputLogFilePath, sep="\n"
  )
}


##################################################################################
##run do_log_reg_chw.R GWAS
##################################################################################
runDoLogRegChwR=function(
  doLogRegChwR=NULL,
  pangenomeFilePath=NULL,
  phenoFilePath=NULL,
  outfix=NULL
  )
{
  gwas = paste0("Rscript ",doLogRegChwR," -bips ",pangenomeFilePath," -phenotype ",phenoFilePath," -prefix ",outfix)
  cat(paste0(gwas,"\n"))
  system(gwas)
  message("logistic regression completed!")
  cat(paste0("Logistic regression GWAS completed: ",Sys.time()), 
      file = inputLogFilePath, append=TRUE, sep="\n")  
}

##################################################################################
##make gemma files
#create gemma_gen_format.txt and gemma_snp_format.txt files for lmm GWAS
#gemma_gen_format = ids, 1, 0, genotypes (1/0)
#gemma_snp_format = snp_id, position, chr (24 for haploids)
##################################################################################
makeGemmaFiles=function(
  pangenome=NULL,
  prefix=NULL
  )
  {
  #gemma_gen
  gemmaGenPath=paste0(prefix,".gemma_gen.txt")
  gen=cbind(row.names(pangenome), rep(1,nrow(pangenome)),rep(0,nrow(pangenome)),pangenome)
  write.table(gen,gemmaGenPath,row.names=F, col.names=F, quote=F)
  #gemma_snp
  gemmaSnpPath=paste0(prefix,".gemma_snp.txt")
  snp=cbind(row.names(pangenome),seq(1:nrow(pangenome)),rep(24, nrow(pangenome)))
  write.table(snp,gemmaSnpPath,row.names=F, col.names=F, quote=F)
  return(c(gemmaGenPath,gemmaSnpPath))
}
##################################################################################
##function for getting log null value from gemma log file
##################################################################################
getLogNull=function(datafiles=NULL){ 
  a=scan(datafiles,what=character(0),sep="\n")
  b=unlist(strsplit(a[17]," "))
  b=as.numeric(b[length(b)])
  return(b)
}

##################################################################################
##function for correcting p-values output from gemma
##################################################################################
getPvalue=function(
  gemma=NULL,
  logNull=NULL)
  {
  D=2*abs(logNull-as.numeric(gemma[length(gemma)])) # This is the LH1 column
  pval=-log10(exp(1))*pchisq(D,1,low=F,log=TRUE)
  return(pval)
}

##################################################################################
##run gemma /dipro/mmm/gorm/v3/ana/Cdif/earle/SNPanalysis/gemma/gemma
##################################################################################
runGemma=function(
  gemmaPath=NULL,
  gemmaFilePaths=NULL,
  phenoFilePath=NULL,
  relm=NULL,
  outfix=NULL
  
  )
{ 
  checkExistence(relm)
  gemmaGenPath=gemmaFilePaths[1]
  gemmaSnpPath=gemmaFilePaths[2]
  if(!file.exists(gemmaGenPath)){
    stop("Gemma gen file is not created properly.")
  }
  if(!file.exists(gemmaSnpPath)){
    stop("Gemma snp file is not created properly.")
  }
lmm_gwas=paste0(gemmaPath," -g ",gemmaGenPath," -p ",phenoFilePath," -a ",gemmaSnpPath," -k ",relm," -lmm 4 "," -o ",paste0(prefix,".lmm")," -maf 0")
cat(paste0(lmm_gwas,"\n"))
system(lmm_gwas)
message("GEMMA completed!")
cat(paste0("GEMMA completed: ",Sys.time()), 
    file = inputLogFilePath, append=TRUE, sep="\n")    
}

##################################################################################
##make manhattan plot
##################################################################################
plotManhattan=function(
  pvalue=NULL,
  prefix=NULL
  )
{
  SNPid = 1:length(pvalue)
  png(paste0(prefix,"_manhattan.png"),width=1200,height=750)
  plot(x = SNPid, y = pvalue, type="p",main=paste0(prefix," GWAS manhattan plot"),xlab="SNPid",ylab="-log10(p)")
  dev.off()
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
pangenome = extractInputArgument(argName="pangenome", commandLineInputs = inputs, default = c(0))
phenotype = extractInputArgument(argName="phenotype", commandLineInputs = inputs, default = c(0))
gemma= extractInputArgument(argName="gemma", commandLineInputs = inputs, default = "no")
relm=extractInputArgument(argName="relm", commandLineInputs = inputs, default = c(0))
prefix = extractInputArgument(argName="prefix", commandLineInputs = inputs, default = c(0))
scriptDir=extractInputArgument(argName="scriptDir", commandLineInputs = inputs, default = getwd())
externalSoftwarePaths = extractInputArgument(argName="externalSoftwarePaths", commandLineInputs = inputs)
runDir= getwd()

#check that gemma argument is correct format
gemma=tolower(gemma)
if(gemma!= "yes" & gemma !="no"){
  stop("\nIncorrect format: -gemma [yes] [no]\n")
}

#check existence of data files and load them in##
checkExistence(pangenome)
pangenome.df=read.table(pangenome, header=T, as.is=T)
checkExistence(phenotype)
phenotype.df = read.table(phenotype, header=T, as.is=T)

#Get the path of external software used ##
externalSoftwarePaths.df = read.table(file=externalSoftwarePaths, header=T, as.is = T) 
gemmaPath = getSoftwarePath("Gemma", externalSoftwarePaths.df)

# Check all called scripts exist #
cat(paste0("Script directory: ",scriptDir))
doLogRegChwR = formatPath(scriptDir,"do_logreg_chw.R")

##### Save input and later run information #####
inputLogFilePath = paste(prefix, "GWAS_logfile.txt", sep="_")
displayInputInfo(
  pangenome = pangenome,
  phenotype = phenotype,
  prefix = prefix,
  gemma = gemma,
  relm = relm,
  scriptDir=scriptDir,
  runDir = runDir,
  inputLogFilePath = inputLogFilePath
)

message(paste(c("Saved input and run information to ",inputLogFilePath, "."), collapse=""))

#order genomes in pan-genome matrix to match order of genomes in relatedness matrix/phenotype file
genomeIds=names(pangenome.df)
if (all(genomeIds!=phenotype.df$name)){
  if(length(genomeIds)==length(phenptype.df$name)){
pangenome.df = pangenome.df[,match(phenotype.df$name, genomeIds)]
}else (stop("Number of isolates in phenotype file does not equal number of isolates in pangenome file!"))
}
#remove genome ids from file (needed for running do_logreg_chw.R)##
bipsFilePath=paste0(pangenome,".bipsFile")
write.table(pangenome.df, bipsFilePath, row.names=T,col.names=F,quote=F,sep="\t")

#run logistic regression and lmm GWAS for given set of phenotypes and pan-genome set##
for (p in 2:ncol(phenotype.df)){
  #create phenotype files
  pheno=phenotype.df[p]
  name=names(pheno)
  outfix = paste0(name,"_",prefix)
  cat(paste0("Doing: ",name,"\n"))
  phenoFilePath = paste0(name,".pheno")
  write.table(pheno,phenoFilePath, row.names=F, col.names=F,quote=F)
  
  #run logistic regression GWAS
  runDoLogRegChwR(doLogRegChwR,bipsFilePath,phenoFilePath,outfix)
  logreg=read.table(paste0(outfix,"_SNP_logreg_output.txt"),header=T,as.is=TRUE,sep="\t")
  plotFix=paste0(outfix,"_SNP_logreg")
  plotManhattan(logreg$Pvalue,plotFix)
  
  #check if gemma option specified and if so run GEMMA##
  if(gemma=="yes"){
    #create gemma_gen_format.txt and gemma_snp_format.txt files for lmm GWAS##
    gemmaFilePaths=makeGemmaFiles(pangenome.df,prefix)
    #run GEMMA
    runGemma(gemmaPath,gemmaFilePaths,phenoFilePath,relm,outfix)

  #correct GEMMA p-values and calculate -logl0 pvalues
  gemmaLogPath=paste0("./output/",prefix,".lmm.log.txt")
  gemmaResPath=paste0("./output/",prefix,".lmm.assoc.txt")
  logNull=getLogNull(gemmaLogPath)
  gemmaResult=read.table(gemmaResPath,header=TRUE,as.is=TRUE)
  negLog10=apply(gemmaResult,1,getPvalue,logNull=logNull)
  gemmaResultCorrected=cbind(gemmaResult,negLog10)
  write.table(gemmaResultCorrected,file=paste0(outfix,"_LMM_corrected.txt"),row=F,col=T,sep="\t",quote=FALSE)
  }
  logreg=read.table(paste0(outfix,"_LMM_corrected.txt"),header=T,as.is=TRUE,sep="\t")
  plotFix=paste0(outfix,"_LMM_corrected")
  plotManhattan(logreg$negLog10,plotFix)
}  
