# File SNP_MAIN_GWAS.R
# Authors: Earle, S., Wu, C.-H. and Wilson, D. J. 
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
	"SNP_GWAS_SVN.R performs a reference based SNP GWAS",
	"Authors: Earle, S., Wilson, D. J. and Wu, C.-H.",
	"Usage: Rscript SNP_GWAS_SVN.R -fastafiles -phenotype -ref_fa -ref_gbk -other",
	"       E.g. Rscript SNP_GWAS_SVN.R -fastafiles Cdif_fasta.txt -phenotype pheno.txt -phylogeny tree.raxml -ref_fa referance.fa etc.",
	"","Options:",
	"-dataFile",
	" File containing list of file paths to mapped fasta files, corresponding ids and phenotypes.",
	"-phylogeny",
	" Input a phylogeny file (Default: builds a phylogeny in PhyML or RAxML)",
	"-ref_fa",
	" Reference fasta file",
	"-ref_gbk",
	" Reference genbank file",
	"-prefix",
	" Output file name prefix",
	"-script_dir",
	" Absolute path of the directory of the location of the scripts (Default: Current working directory)",
	"-CFML_prefix",
	" If ClonalFrameML has previously been run, the prefix of the file names (Default: run ClonalFrameML)",
	"-run_gemma",
	" Option to additionally run GEMMA (yes or no) Default: no",
	"-PCA",
	" Option to additionally run a PCA analysis using Eigensoft (yes or no)",
	"-n_PCs",
	" If running PCA, how many principle components to compute (Default: 100)",
	"-externalSoftwarePaths",
	"A comma delimited file containing the name and paths of the external software used in the analysis",
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
## Check the file path.
##################################################################################
# checkExistence = function(objectName = NULL, path = NULL){
	# if(!file.exists(path)){
		# if(is.null(objectName)){
			# stop(paste(c("The file",path,"does not exist!"), collapse=" ",sep=""))
		# }else{
			# stop(paste(c("The file", objectName, "is not found at", path ,"!"), collapse=" ",sep=""))			
		# }
	# }
	
# }

##################################################################################
## Turn command line inputs into a n x 2 matrix
##################################################################################
getCommandLineInputMatrix = function(args = NULL){
	argsCount = length(args)/2
	inputs = matrix(nrow=argsCount, ncol=2)
	inputs[,1] = args[c(1:argsCount)*2 -1]
	
	if(length(which(regexpr("-",inputs[,1]) != 1)) > 0){
		for(i in 1:nrow(inputs)){
			message(paste(c(i,": ", inputs[i,1],", ", regexpr("-",inputs[i,1]) ), collapse="") )
		}
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
## Varify the existence of a file.
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
	dataFile = NULL, 
	phylogeny = NULL, 
	prefix = NULL, 
	ref_fa = NULL,
	ref_gbk = NULL,
	script_dir = NULL,
	inputLogFilePath = NULL){
	if(!is.null(phylogeny)){
		cat(paste0("Data file: ", dataFile),
			paste0("Output file name prefix: ",prefix),
			paste0("Phylogeny: ",phylogeny),
			paste0("Reference fasta file: ",ref_fa),
			paste0("Reference genbank file: ",ref_gbk),
			paste0("Script directory: ",script_dir),
			paste0("Start time: ", Sys.time()),
			file = inputLogFilePath, sep="\n"
		)
	}else{
		cat(paste0("Data file: ", dataFile),
			paste0("Output file name prefix: ",prefix),
			paste0("Phylogeny: Will build phylogeny in either RAxML or PhyML"),
			paste0("Script directory: ",script_dir),
			paste0("Start time: ",Sys.time()),
			file = inputLogFilePath, sep="\n"
		)
	}
}

##################################################################################
## Get prefix for CFML files.
##################################################################################
getCfmlPrefix = function(cfmlPrefix = NULL, prefix = NULL){	
	if(cfmlPrefix!=0){	
		checkCfmlInputs(cfmlPrefix)
		return(cfmlPrefix)		
	}else{
		return(prefix)		
	}	
}

##################################################################################
## Check the existence of CFML files.
##################################################################################
checkCfmlInputs = function(cfmlPrefix = NULL){	
	
	cfmlpath=unlist(strsplit(cfmlPrefix,"/"))
	cfmlname=cfmlpath[length(cfmlpath)]
	cfmlpath=paste(cfmlpath[1:(length(cfmlpath)-1)],collapse="/")

	if((cfmlpath==cfml_prefix)==TRUE){
		existance =dir("./",pattern=glob2rx(paste0(cfml_prefix,"*")))
		
		if(length(existance)<2){
			stop("\nIncorrect usage: CFML files do not exist\n")
		}
		
	} else {
		existance =dir(cfmlpath,pattern=glob2rx(paste0(cfmlname,"*")))

		if(length(existance)<2){
			stop("\nIncorrect usage: CFML files do not exist\n")
		}
	}
	
}

##################################################################################
## Get the phylogeny for creating CFML output files.
##################################################################################
getPhylogenyForCFML = function(
  phylogeny = NULL,
  fastafiles = NULL,
  raxmlBuildR = NULL,
  raxmlPath = NULL,
  phymlBuildR = NULL,
  phymlPath = NULL,
  ref_fa = NULL,
  dataFile = NULL,
  prefix = NULL){
  
  ## If no phylogeny inputted, build tree.
  if(!is.null(phylogeny)){
    phylogeny = phylogeny
    
    ## Check existance of tree file
    existance = file.exists(phylogeny)
    if(length(which(existance==FALSE))>0) {
      stop("\nIncorrect usage: phylogeny file doesn't exist\n")
    }
    
  } else if(is.null(phylogeny) & length(fastafiles)>=100){
  	message("Running RAxML")
    
    phylogeny = buildTreeByRAxML(raxmlBuildR = raxmlBuildR, raxmlPath = raxmlPath, ref_fa = ref_fa, dataFile = dataFile, prefix = prefix)
    
  } else if(is.null(phylogeny) & length(fastafiles)<100){
    message("Running PhyML")
    phylogeny = buildTreeByPhyML(phymlBuildR = phymlBuildR, phymlPath = phymlPath, ref_fa = ref_fa, dataFile = dataFile, prefix = prefix)
    
  }
  
  return(phylogeny)
  
}

##################################################################################
## Build the tree using RAxML.
##################################################################################
buildTreeByRAxML = function(raxmlBuildR = NULL, 
                            raxmlPath = NULL,
                            ref_fa = NULL, 
                            dataFile = NULL, 
                            prefix = NULL){
  
  cat(paste0("Building tree using RAxML"), sep="\n")
  system(paste(c("Rscript",raxmlBuildR, ref_fa, dataFile, prefix, raxmlPath, FALSE), collapse=" "))
  phylogeny = paste0("RAxML_bestTree.",prefix)
  cat(paste0("Phylogeny built in RAxML completed: ", Sys.time()),
      file=paste0(prefix,"_input_logfile.txt"), append=TRUE, sep="\n")
  raxmlInfoPath = paste("RAxML_info",prefix,sep=".")
  
  return(phylogeny)
  
}

##################################################################################
## Build the tree using PhyML.
##################################################################################
buildTreeByPhyML = function(phymlBuildR = NULL,
                            phymlPath = NULL,
                            ref_fa = NULL, 
                            dataFile = NULL, 
                            prefix = NULL){
  
  cat(paste0("Building tree using PhyML"),sep="\n")
  system(paste(c("Rscript", phymlBuildR, ref_fa, dataFile, prefix, phymlPath, FALSE), collapse=" "))
  phylogeny=paste(prefix,"_phylip.phylip_phyml_tree.txt",sep="")
  cat(paste0("Phylogeny built in PhyML completed: ",Sys.time()),
      file=paste0(prefix,"_input_logfile.txt"),append=TRUE,sep="\n")
  
  return(phylogeny)
  
}


##################################################################################
## Run snp_calling.R.
##################################################################################
runSnpCalling = function(snpCallingR = NULL,
	prefix = NULL,
	ref_fa = NULL,
	dataFile = NULL,
	inputLogFilePath = NULL){
		
	system(paste(c("Rscript", snpCallingR, prefix, ref_fa, dataFile), collapse = " "))
	cat(paste0("SNP calling completed: ",Sys.time()), file=inputLogFilePath, append=TRUE, sep="\n")
	
	bipinfo = paste(c(formatDir(getwd()), prefix, ".bipinfo.txt"), collapse="")
	bippatterns = paste(c(formatDir(getwd()), prefix, ".bip.patterns.txt"), collapse="")
	snpinfo = paste(c(formatDir(getwd()), prefix, ".snpinfo.txt"), collapse="")
	snppatterns = paste(c(formatDir(getwd()), prefix, ".snp.patterns.txt"), collapse="")
	checkExistence(filePaths=bipinfo, msg="The bipinfo file is not properly created by snp_calling.R")
	checkExistence(filePaths=bippatterns, msg="The bippatterns file is not properly created by snp_calling.R")
	checkExistence(filePaths=snpinfo, msg="The snpinfo file is not properly created by snp_calling.R")
	checkExistence(filePaths=snppatterns, msg="The snppatterns file is not properly created by snp_calling.R")
	return(list(bipinfo = bipinfo, bippatterns = bippatterns, snpinfo = snpinfo, snppatterns = snppatterns))

}

##################################################################################
## Run get_annotation.R.
##################################################################################
runGetAnnotation = function(
	getAnnotationR = NULL,
	bipinfo = NULL,
	ref_fa = NULL,
	ref_gbk = NULL,
	prefix = NULL,
	inputLogFilePath = NULL){
		
	system(paste(c("Rscript", getAnnotationR, bipinfo, ref_fa, ref_gbk, prefix), collapse = " "))
	cat(paste0("Biallelic SNP annotation completed: ",Sys.time()), file = inputLogFilePath, append = TRUE, sep="\n")

	bipAnnotationFile = paste(c(formatDir(getwd()), prefix, ".bip.annotation.txt"), collapse="")
	checkExistence(bipAnnotationFile)
	return(bipAnnotationFile)
	
}

##################################################################################
## Generate input files for imputation from CFML.
##################################################################################
generateCfmlInputs = function(createFastaPath = NULL,
                              prefix = NULL,
                              runCfmlR = NULL,
                              clonalFrameMLPath = NULL,
                              phylogeny = NULL,
                              cfmlFastaFilePath = NULL,
                              dataFile = NULL,
                              getMissingSiteCountR = NULL){
  
  
  cfmlFastaFilePath = runCreateFasta(
      createFastaPath = createFastaPath, 
      prefix = prefix,
      dataFile = dataFile
  )
  message("phylip to fasta version completed!")
  message(paste(c("phylogeny:", phylogeny), collapse=" "))
 
  cfmlOutputPaths = runCFML(runCfmlR = runCfmlR,
                            clonalFrameMLPath = clonalFrameMLPath,
                            phylogeny = phylogeny,
                            cfmlFastaFilePath = cfmlFastaFilePath,
                            dataFile = dataFile,
                            getMissingSiteCountR = getMissingSiteCountR,
                            prefix = prefix)
  

  
  return(list(cfmlPositionCrossRef = cfmlOutputPaths$positionCrossRef,
              cfmlMLSeqFasta = cfmlOutputPaths$mlSeqFasta))
  
}


##################################################################################
## Get ClonalFrameML inputs.
##################################################################################
getCfmlInputs = function(cfml_prefix = NULL,
                         prefix = NULL,
                         fastafiles = NULL,
                         raxmlBuildR = NULL, 
                         raxmlPath = NULL,
                         phymlBuildR = NULL,
                         phymlPath = NULL, 
                         ref_fa = NULL, 
                         dataFile = NULL,
                         createFastaPath = NULL,
                         runCfmlR = NULL,
                         clonalFrameMLPath = NULL,
                         phylogeny = NULL,
                         cfmlFastaFilePath = NULL,
                         getMissingSiteCountR = NULL){
  
  if(cfml_prefix==prefix){
    phylogeny = getPhylogenyForCFML(phylogeny = phylogeny, fastafiles = fastafiles, raxmlBuildR = raxmlBuildR, raxmlPath = raxmlPath,
                        phymlBuildR = phymlBuildR, phymlPath = phymlPath, ref_fa = ref_fa, dataFile = dataFile, prefix = prefix)
    
    message("phyML completed!")
    message("createFastaPath: ", createFastaPath)
    
    cfmlInputPaths = generateCfmlInputs(createFastaPath = createFastaPath, prefix = prefix, runCfmlR = runCfmlR,
                                        clonalFrameMLPath = clonalFrameMLPath, phylogeny = phylogeny, cfmlFastaFilePath = cfmlFastaFilePath,
                                        dataFile = dataFile, getMissingSiteCountR = getMissingSiteCountR)
    
    message("clonalFrameML completed!")
    
  }else{
    cfmlInputPaths = list(cfmlPositionCrossRef = paste(cfml_prefix,".position_cross_reference.txt",sep=""),
                          cfmlMLSeqFasta = paste(cfml_prefix,".ML_sequence.fasta",sep=""))    
  }
  return(cfmlInputPaths)
}

##################################################################################
## Run ConvertPhylipToFasta.
##################################################################################
runCreateFasta = function(createFastaPath = NULL, prefix = NULL, dataFile = NULL){
		
	cfmlFastaFilePath =  paste(prefix,".fasta", sep="")
	pathSplit = unlist(strsplit(createFastaPath, split="/"))
	classPath = gsub(pathSplit[length(pathSplit)], "", createFastaPath)
	createFasta = gsub(".class", "", pathSplit[length(pathSplit)])
	createFastaCommand = paste(c("java", "-cp", classPath, createFasta, prefix, dataFile), collapse=" ")
	
	system(createFastaCommand)
	checkExistence(cfmlFastaFilePath, "The fasta file has not been correctly created by ConvertPhylipToFasta!")
	return(cfmlFastaFilePath)
	
}

##################################################################################
## Run runCFML.R.
##################################################################################
runCFML = function(runCfmlR = NULL,
	clonalFrameMLPath = NULL,
	phylogeny = NULL,
	cfmlFastaFilePath = NULL,
	kappa = 2.0,
	dataFile = NULL,
	getMissingSiteCountR = NULL,
	prefix = NULL){
  
    message(paste(c("runCfmlR:",runCfmlR), collapse=" "))
		runCFMLCommand = paste(c(
			"Rscript", runCfmlR, 
			clonalFrameMLPath, 
			phylogeny, 
			cfmlFastaFilePath,
			2.0, #default kappa value in ClonalFrameML
			dataFile,
			getMissingSiteCountR,
			prefix),
		collapse=" ")
	system(runCFMLCommand)
	
	

	positionCrossRef = paste(c(formatDir(getwd()), prefix,".position_cross_reference.txt"), collapse="")
	mlSeqFasta = paste(c(formatDir(getwd()), prefix,".ML_sequence.fasta"), collapse="")
	
	checkExistence(positionCrossRef)
	checkExistence(mlSeqFasta)		
	return(list(positionCrossRef = positionCrossRef, mlSeqFasta = mlSeqFasta))
	
}

##################################################################################
## Run Regression_format.R.
##################################################################################
runRegressionFormat = function(
	prefix = NULL, 
	regressionFormatR = NULL, 
	bipinfo = NULL,
	bippatterns = NULL,
	positionCrossRef = NULL,
	mlSeqFasta = NULL,
	dataFile = NULL){
		
	#positionCrossRef = paste(prefix,".position_cross_reference.txt", sep="")
	#mlSeqFasta = paste(prefix,".ML_sequence.fasta", sep="")
	system(paste(c("Rscript", regressionFormatR, positionCrossRef, mlSeqFasta, bipinfo, bippatterns, prefix, dataFile), collapse=" "))

	gemmaGenPath = paste(c(formatDir(getwd()), prefix, "_gemma_gen_format.txt"), collapse="")
	gemmaSnpPath = paste(c(formatDir(getwd()), prefix, "_gemma_snp_format.txt"), collapse="")
	if(!file.exists(gemmaGenPath)){
		stop("Gemma gen file is not created properly.")
	}
	if(!file.exists(gemmaSnpPath)){
		stop("Gemma snp file is not created properly.")
	}
	return(list(gemmaGenPath = gemmaGenPath, gemmaSnpPath = gemmaSnpPath))
	
}

##################################################################################
## Run snpLogRegression.R.
##################################################################################
runSnpLogRegression = function(gemmaGenPath = NULL,
	gemmaSnpPath = NULL,
	prefix = NULL,
	dataFile = NULL,
	annotationPath = NULL,
	gwasPlotsR = NULL,
	inputLogFilePath = NULL){
		
	system(paste(c("Rscript", snpLogRegressionR, 
				gemmaGenPath, gemmaSnpPath, prefix, dataFile, annotationPath, gwasPlotsR), 
				collapse=" "))
				
	cat(paste0("Logistic regression for biallelic SNPs completed: ",Sys.time()), 
		file = inputLogFilePath, append=TRUE, sep="\n")

}

##################################################################################
## Run tri_tetra_allelic_logreg.R.
##################################################################################
runTriTetraAllelicLogreg = function(triTetraAllelicLogregR = NULL, 
	cfmlPositionCrossRef = NULL, 
	cfmlMLSeqFasta = NULL, 
	prefix = NULL, 
	dataFile = NULL, 
	snpinfo = NULL, 
	snppatterns = NULL,
	annotFilePath = NULL,
	inputLogFilePath = NULL){
	message("flag1")	
	message(cfmlPositionCrossRef)
	message(cfmlMLSeqFasta)
	triTetraAllelicLogregCommand = paste(c("Rscript", triTetraAllelicLogregR, 
			cfmlPositionCrossRef, 
			cfmlMLSeqFasta, 
			prefix, 
			dataFile, 
			snpinfo, 
			snppatterns,
			annotFilePath), 
		collapse=" "
	)
	
	message("flag2")
	message(triTetraAllelicLogregCommand)
	system(triTetraAllelicLogregCommand)
	message("flag3")
	cat(paste0("Logistic regression for tri and tetra alellic SNPs completed: ",Sys.time()),
		file=paste0(prefix,"_input_logfile.txt"), append=TRUE,sep="\n"
	)
	
}

##################################################################################
## Run get_annotations_tritetra_chw.R
##################################################################################
runGetAnnotationsTriTetra = function(
	getAnnotationsTriTetra = NULL,
	snpinfo = NULL,
	ref_fa = NULL,
	ref_gbk = NULL,
	prefix = NULL,
	inputLogFilePath = NULL){
		
	system(paste(c("Rscript", getAnnotationsTritetraR, snpinfo, ref_fa, ref_gbk, prefix), collapse=" "))
	snpAnnotPath = paste(c(formatDir(getwd()),prefix,".snp.annotation.txt"),collapse="")
	message(snpAnnotPath)
	checkExistence(snpAnnotPath, "The snp annotation file has not been created correctly by the get_annotations_tritetra_chw.R")
	cat(paste0("Tri and tetra allelic SNP annotation completed: ",Sys.time()),
		file=inputLogFilePath, append=TRUE, sep="\n")
	return(snpAnnotPath)
		
	
}

##################################################################################
## Run snpGemma.R
##################################################################################
runSnpGemma = function(snpGemmaR = NULL,
	prefix = NULL,
	gemmaPath = NULL,
	gemmaGenPath = NULL,
	gemmaSnpPath = NULL,
	dataFile = NULL,
	bipAnnotPath = NULL,
	signif.threshold = NULL,
	gwasPlotsR = NULL,
	inputLogFilePath = NULL){
	message("Run snpGemma")

				
	snpGemmaCommand = paste(c("Rscript", snpGemmaR,
		prefix, gemmaPath, gemmaGenPath, gemmaSnpPath, dataFile, bipAnnotPath, signif.threshold, gwasPlotsR), collapse=" ")
	system(snpGemmaCommand)
	cat(paste0("GEMMA completed: ",Sys.time()), file = inputLogFilePath, append=TRUE, sep="\n")
	
}

##################################################################################
## Run EigenSoft.R
##################################################################################
runEigenSoft = function(eigenSoftFormatR = NULL,
	eigenSoftPath = NULL,
	gemmaGenPath = NULL,
	gemmaSnpPath = NULL,
	prefix = NULL,
	dataFile = NULL,
	npc = NULL){
	
	eigenSoftFormatCommand = paste(c("Rscript", eigenSoftFormatR, 
		eigenSoftPath, gemmaFilePaths$gemmaGenPath, gemmaFilePaths$gemmaSnpPath, prefix, dataFile, npc), 
		collapse=" ")
	system(eigenSoftFormatCommand)
	evecFile = paste(prefix,"pca.evec", sep=".")
	evalFile = paste(prefix,"pca.eval", sep=".")
	checkExistence(evecFile, "The PCA evec file has not been created correctly by EigenSoft.R")	
	cat(paste0("Saved input and run information to ",prefix,"_input_logfile.txt"),sep="\n")
	return(list(evecFile = evecFile, evalFile = evalFile))
}

##################################################################################
## Run snpLogRegressionPCA.R
##################################################################################
runSnpLogRegPca = function(snpLogRegressionPcaR = NULL, 
		gemmaGenPath = NULL, 
		gemmaSnpPath = NULL, 
		eigenSoftScores = NULL,
		eigenSoftValues = NULL,
		prefix = NULL,
		dataFile = NULL,
		bipAnnotationFile = NULL,
		gwasPlotsR = NULL){
  
	snpLogRegPcaCommand = paste(c("Rscript", snpLogRegressionPcaR, 
		gemmaGenPath, 
		gemmaSnpPath, 
		eigenSoftScores,
		eigenSoftValues,
		prefix,
		dataFile,
		bipAnnotationFile,
		gwasPlotsR),
		collapse=" "
	)
	system(snpLogRegPcaCommand)
	
	
}



## Read the command line arguments
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

message("Start reading command inputs ...")
inputs = getCommandLineInputMatrix(args = args)
dataFile = extractInputArgument(argName="dataFile", commandLineInputs = inputs, default = c(0))
phylogeny = extractInputArgument(argName="phylogeny", commandLineInputs = inputs)
ref_fa= extractInputArgument(argName="ref_fa", commandLineInputs = inputs, default = c(0))
ref_gbk= extractInputArgument(argName="ref_gbk", commandLineInputs = inputs, default = c(0))
script_dir= extractInputArgument(argName="script_dir", commandLineInputs = inputs, default = getwd())
cfml_prefix = extractInputArgument(argName="CFML_prefix", commandLineInputs = inputs, default = c(0))
prefix = extractInputArgument(argName="prefix", commandLineInputs = inputs, default = c(0))
gemma = extractInputArgument(argName="run_gemma", commandLineInputs = inputs, default = "no") 
eigensoft = extractInputArgument(argName="PCA", commandLineInputs = inputs, default = "no")
npc = extractInputArgument(argName="n_PCs", commandLineInputs = inputs, default = 100)
externalSoftwarePaths = extractInputArgument(argName="externalSoftware", commandLineInputs = inputs)


## Check that the provided paths to the ClonalFrameML outputs are valid
cfml_prefix = getCfmlPrefix(cfml_prefix, prefix)

## Read in the paths of the fasta files and check that the paths are valid.
checkExistence(dataFile)
data.df = read.table(file = dataFile, header = T, as.is = T)
fastafiles = data.df$filePath
phenotype = data.df$phenotype

## If a phylogenetic has been provided, check whether the path is valid.
if(!is.null(phylogeny)){
	checkExistence(phylogeny)
}

## Check that the fasta file of the reference sequence exists.
checkExistence(ref_fa)

## Check that the genbank file of the reference sequence exists.
checkExistence(ref_gbk)

## Check the input argument of the gemma option.
gemma=tolower(gemma)
if(gemma!= "yes" & gemma !="no"){
	stop("\nIncorrect format: -gemma [yes] [no]\n")
}


##### Check all called scripts exist #####
cat(paste0("Script directory: ",script_dir))
snpCallingR = formatPath(script_dir,"snpCalling.R")
raxmlBuildR = formatPath(script_dir,"RAxML_build.R")
phymlBuildR = formatPath(script_dir,"PhyML_build.R")
regressionFormatR = formatPath(script_dir,"RegressionFormat.R")
getAnnotationR = formatPath(script_dir,"getAnnotations.R")
triTetraAllelicLogregR = formatPath(script_dir,"ttpLogReg.R")
getAnnotationsTritetraR = formatPath(script_dir,"getAnnotationsTTP.R")
eigenSoftFormatR = formatPath(script_dir,"EigensoftFormat.R")
runCfmlR = formatPath(script_dir,"runCFML.R")
getMissingSiteCountR = formatPath(script_dir,"getMissingSiteCount.R")
snpLogRegressionR = formatPath(script_dir,"snpLogRegression.R")
gwasPlotsR = formatPath(script_dir,"gwasPlots.R")
snpGemmaR = formatPath(script_dir,"snpGEMMA.R")
snpLogRegressionPcaR = formatPath(script_dir, "snpLogRegressionPCA.R")


##### Get the path of external software used #####
externalSoftwarePaths.df = read.table(file=externalSoftwarePaths, header=T, as.is = T) 
clonalFrameMLPath = getSoftwarePath("ClonalFrameML", externalSoftwarePaths.df)
gemmaPath = getSoftwarePath("GEMMA", externalSoftwarePaths.df)
raxmlPath = getSoftwarePath("RAxML", externalSoftwarePaths.df)
phymlPath = getSoftwarePath("PhyML", externalSoftwarePaths.df)
createFastaPath = getSoftwarePath("CreateFasta", externalSoftwarePaths.df)
eigenSoftPath = getSoftwarePath("EigenSoft", externalSoftwarePaths.df)

##### Save input and later run information #####
inputLogFilePath = paste(prefix, "input_logfile.txt", sep="_")
displayInputInfo(
	dataFile = dataFile, 
	phylogeny = phylogeny, 
	prefix = prefix, 
	ref_fa = ref_fa,
	ref_gbk = ref_gbk,
	script_dir = script_dir,
	inputLogFilePath = inputLogFilePath
)
message(paste(c("Saved input and run information to ",inputLogFilePath, "."), collapse=""))


##### SNP calling #####
snpCallingOutputPaths = runSnpCalling(snpCallingR = snpCallingR, prefix = prefix, ref_fa = ref_fa, dataFile = dataFile, inputLogFilePath = inputLogFilePath)
bipinfo=snpCallingOutputPaths$bipinfo
bippatterns=snpCallingOutputPaths$bippatterns
snpinfo=snpCallingOutputPaths$snpinfo
snppatterns=snpCallingOutputPaths$snppatterns

cfmlInputNames = getCfmlInputs(cfml_prefix = cfml_prefix, prefix = prefix, fastafiles = fastafiles, 
              raxmlBuildR = raxmlBuildR, raxmlPath = raxmlPath,
              phymlBuildR = phymlBuildR, phymlPath = phymlPath, ref_fa = ref_fa, dataFile = dataFile, 
              createFastaPath = createFastaPath, runCfmlR = runCfmlR, clonalFrameMLPath = clonalFrameMLPath, 
              phylogeny = phylogeny, cfmlFastaFilePath = cfmlFastaFilePath, getMissingSiteCountR = getMissingSiteCountR)
cfmlPositionCrossRef = cfmlInputNames$cfmlPositionCrossRef
cfmlMLSeqFasta = cfmlInputNames$cfmlMLSeqFasta


##### Annotate SNPs with gene names etc. #####
bipAnnotationFile = runGetAnnotation(getAnnotationR, bipinfo, ref_fa, ref_gbk, prefix, inputLogFilePath)

##### Formatting for log reg and gemma and running log reg #####
gemmaFilePaths = runRegressionFormat(prefix, regressionFormatR, bipinfo, bippatterns, cfmlPositionCrossRef, cfmlMLSeqFasta, dataFile)
gemmaGenPath = gemmaFilePaths$gemmaGenPath
gemmaSnpPath = gemmaFilePaths$gemmaSnpPath

##Do raw logisitc regression on binary snp data
runSnpLogRegression(gemmaGenPath, gemmaSnpPath, prefix, dataFile, bipAnnotationFile, gwasPlotsR, inputLogFilePath)

## Annotate SNPs with gene names etc.
annotFilePath = runGetAnnotationsTriTetra(getAnnotationsTritetraR, snpinfo, ref_fa, ref_gbk, prefix, inputLogFilePath)

## Formatting and running tri and tetra allelic site logistic regression
runTriTetraAllelicLogreg(triTetraAllelicLogregR,
	cfmlPositionCrossRef, cfmlMLSeqFasta, prefix, dataFile, snpinfo, snppatterns, annotFilePath, inputLogFilePath)

## Run gemma
if(gemma=="yes"){	
	runSnpGemma(snpGemmaR, 
		prefix, gemmaPath, gemmaGenPath, gemmaSnpPath, dataFile, bipAnnotationFile, 0, gwasPlotsR, inputLogFilePath)
}

##If Eigensoft is yes, run
if(eigensoft=="yes"){
	#Format for eigensoft and run
	eigenFiles = runEigenSoft(eigenSoftFormatR, eigenSoftPath, gemmaGenPath, gemmaSnpPath, prefix, dataFile, npc)
  message("Eigenfiles created.")
	runSnpLogRegPca(snpLogRegressionPcaR, 
		gemmaGenPath, gemmaSnpPath, eigenFiles$evecFile, eigenFiles$evalFile, prefix, dataFile, bipAnnotationFile, gwasPlotsR)
		
}


