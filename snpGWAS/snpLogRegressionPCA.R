###############################################################################
## Do logistic regression
###############################################################################
logreg<-function(
	data = NULL, 
	pheno = NULL, 
	additionalPredictors = NULL){
		
		fullData = data.frame(snp = as.factor(data), additionalPredictors)
		logr1=glm((pheno==1)~., family=binomial(logit),data = additionalPredictors)$deviance
		logr2=glm((pheno==1)~., family=binomial(logit), data = fullData)$deviance
#     message(paste(c("# predictors:", length(additionalPredictors)) , collapse=" "))
#     message(paste(c("logr1:",logr1, collapse=" ")))
# 		message(paste(c("logr2:",logr2, collapse=" ")))
 		dev=abs(diff(c(logr1,logr2)))
		logchi = -log10(pchisq(dev,1,low=F))
		return(logchi)
}


###############################################################################
## Calculate the odds from a 2 x 2 contingency table
###############################################################################
calculateOdds<-function(tb = NULL){
	#control = tb[1,]
	#cases = tb[2,]
	if(!(any(tb['0',]==0) && any(tb['1',]==0))){
	#if(length(which(tb[1,]==0)) ==0 | length(which(tb[2,]==0)) ==0){
		odds=(tb[2,2]*tb[1,1])/(tb[1,2]*tb[2,1])
	} else {
		odds="NA"
	}
	
}


###############################################################################
## Perform logistic regression, calculate odds ratio and produce summary for
## a give snp.
###############################################################################
logRegOutput<-function(
	snpData = NULL,
	phenotype = NULL,	
	additionalPredictors = NULL,
	logreg_output_filename = NULL,
	snpPosition = NULL,
	iteration = NULL){
	
	logregression = logreg(
		data = snpData,
		pheno = phenotype,
		additionalPredictors = additionalPredictors
	)
	
	tb = table(phenotype, snpData)
	odds = calculateOdds(tb)


	createLogRegOutputFile(
		iteration = iteration,
		snpPosition = snpPosition, 
		contingency.tb = tb, 
		odds = odds,
		neglog10P = logregression,
		logreg_output = logreg_output_filename
	)

}

###############################################################################
## Create output which summarises the statistical analysis
###############################################################################
createLogRegOutputFile<-function(
	iteration = NULL,
	snpPosition = NULL,
	contingency.tb = NULL,
	odds = NULL,
	neglog10P = NULL,	
	logreg_output_filename = NULL){
	
	info3 = c(paste("SNP", iteration, sep=""), snpPosition)	
	
	logregSummary<-c(info3,
			contingency.tb['0',], #control
			contingency.tb['1',], #cases
			odds,
			neglog10P)
	#message("Hello!")		
	cat(paste(logregSummary, collapse="\t"),
		file = logreg_output_filename, 
		append = T, 
		sep="\n"
	)
	#message("bye")
}

getNumberOfPCs = function(evalFile = NULL){
  eval.df = read.table(file=evalFile, header=F)
  clustering = kmeans(eval.df[,1],2)
  index = which(clustering$cluster==clustering$cluster[1])
  return(index)
}




snpAnalysisLogReg<-function(
	gemmaGenInputPath = NULL,
	gemmaSNPInputPath = NULL,
	eigenSoftScoresPath = NULL,
	eigenSoftValuesPath = NULL,
	prefix = NULL,
	filepathsPhenoTablePath = NULL){
		
	message("snpAnalysisLogReg:")
	message(eigenSoftScoresPath)
	message(eigenSoftValuesPath)
	pcIndex = getNumberOfPCs(eigenSoftValuesPath)
  message(paste(c("#Evec:", length(pcIndex)), collapse=" "))
  
	additionalPredictor.df<-read.table(file=eigenSoftScoresPath, header=F)
  pcCount = (ncol(additionalPredictor.df) - 2)
  if(pcCount < length(pcIndex)){
    pcIndex = c(1:pcCount)    
  }
  message(paste(pcIndex))
	
	message("flag2")

  if(length(pcIndex)==1){
    additionalPredictor.df<-data.frame(additionalPredictor.df[,pcIndex+1])
  }else{
    additionalPredictor.df<-additionalPredictor.df[,pcIndex+1]
  }
	
	colnames(additionalPredictor.df) = paste("pca", c(pcIndex), sep="")
	
	filepathsPhenoTable<-read.table(file=filepathsPhenoTablePath, header=T)
	#message(names(filepathsPhenoTable))
		
	#Get file paths	
	filepaths=filepathsPhenoTable$filePath
	#Get phenotype
	phenotype=as.numeric(filepathsPhenoTable$phenotype)
	

	logreg_output_filename = paste(prefix,"_logreg_pca_output.txt",sep="")
	
	gemmaGenInputConn<-file(gemmaGenInputPath, open = "r")
	gemmaSNPInputConn<-file(gemmaSNPInputPath, open = "r")
	
	cat(paste(c("snpID","position", "control0", "control1", "case0", "case1", "oddsRatio", "negLog10P"), collapse="\t"),
		file = logreg_output_filename, 
		append = F, 
		sep="\n"
	)
	
	#snp.count = nrow(bipinfo)
	#message(paste("Total number of binary snps: ",snp.count))
	iteration = 0
	while (length(gemmaGenLine <- readLines(gemmaGenInputConn, n = 1, warn = FALSE)) > 0) {
		iteration = iteration + 1
		gemmaSNPLine <- readLines(gemmaSNPInputConn, n = 1, warn = FALSE)
		
		#Read binary data
		snpData<-unlist(strsplit(gemmaGenLine, "\t"))
		snpData<-as.numeric(snpData[4:length(snpData)])
		
		## get bip position
		snpPosition = as.numeric(unlist(strsplit(gemmaSNPLine, "\t"))[2])	
		
		
		logRegOutput(
			snpData = snpData,
			phenotype = phenotype,	
			additionalPredictors = additionalPredictor.df,
			logreg_output_filename = logreg_output_filename,
			snpPosition = snpPosition,
			iteration = iteration
		)
	
			
		if((iteration%%1000)==0) cat("Done", iteration ,"\n")
	}
	
	logreg_output_path = paste(getwd(),"/", logreg_output_filename,sep="",collapse="")
	return(logreg_output_path) 
}


runLogReg<-function(
	gemmaGenInputPath = NULL,
	gemmaSNPInputPath = NULL,
	eigenSoftScoresPath = NULL,
	eigenSoftValuesPath = NULL,
	prefix = NULL,
	filepathsPhenoTablePath = NULL,
	annotationPath = NULL){
  message("runLogReg")
			
	logRegOutputPath<-snpAnalysisLogReg(
		gemmaGenInputPath = gemmaGenInputPath,
		gemmaSNPInputPath = gemmaSNPInputPath,
		eigenSoftScoresPath = eigenSoftScoresPath,
		eigenSoftValuesPath = eigenSoftValuesPath,
		prefix = prefix,
		filepathsPhenoTablePath = filepathsPhenoTablePath
	)



	# Read logistic regression output
	logreg=read.table(logRegOutputPath,header=T, as.is=T)

	position = logreg$position
	pvalue = as.numeric(logreg$negLog10P)
	
	manhattan.plot(
		position = position, 
		pvalue = pvalue,
		file = paste(prefix,"_SNP_manhattan.png",sep=""),
		plot.title = paste(prefix," SNP GWAS logistic regression manhatten plot",sep="")
	)	
	
	annot=read.table(file = annotationPath, header=T, sep="\t", quote="")	
	logregSNPAnalysisInfo<-cbind(logreg, annot)
	
	pValOrder = order(logreg$negLog10P,decreasing=T)
	write.table(logregSNPAnalysisInfo[pValOrder,], file=paste(prefix,"_logreg_pca_annot.txt",sep=""), quote=F, row=F, sep="\t")
	
}

args = commandArgs(trailingOnly = TRUE)

if(!(length(args) == 8 | length(args) ==0)) {
  
		message(length(args))
		message(length(args) == 8)
        stop("\nIncorrect usage\n")
}

source(args[8])

runLogReg(
	gemmaGenInputPath = args[1],
	gemmaSNPInputPath = args[2],
	eigenSoftScoresPath = args[3],
	eigenSoftValuesPath = args[4],
	prefix = args[5],
	filepathsPhenoTablePath = args[6],
	annotationPath = args[7]
)

