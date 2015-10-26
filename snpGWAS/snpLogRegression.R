
###############################################################################
## Do logistic regression
###############################################################################
logreg=function(data,pheno){
		logr1=glm((pheno==1)~1,family=binomial(logit))$deviance
		logr2=glm((pheno==1)~1+data,family=binomial(logit))$deviance
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
	logreg_output_filename = NULL,
	snpPosition = NULL,
	iteration = NULL){
		
		
	
	logregression=logreg(snpData,phenotype)

	tb = table(phenotype, snpData)
	odds = calculateOdds(tb)

	createLogRegOutputFile(
		iteration = iteration,
		snpPosition = snpPosition, 
		contingency.tb = tb, 
		neglog10P = logregression, 
		odds = odds,
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
	neglog10P = NULL,
	odds = NULL,
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
		append = TRUE, 
		sep="\n"
	)
	#message("bye")
}



snpAnalysisLogReg<-function(
	gemmaGenInputPath = NULL,
	gemmaSNPInputPath = NULL,
	prefix = NULL,
	filepathsPhenoTablePath = NULL){
		
	filepathsPhenoTable<-read.table(file=filepathsPhenoTablePath, header=T)
	#message(names(filepathsPhenoTable))
		
	#Get file paths	
	filepaths=filepathsPhenoTable$filePath
	#Get phenotype
	phenotype=as.numeric(filepathsPhenoTable$phenotype)
	

	logreg_output_filename = paste(prefix,"_logreg_output.txt",sep="")
	cat(paste(c("id","position","ctrl0","ctrl1","case0","case1","oddsRatio","negLog10PValue"), collapse="\t"),
		file = logreg_output_filename, append=FALSE, sep="\n")
	
	gemmaGenInputConn<-file(gemmaGenInputPath, open = "r")
	gemmaSNPInputConn<-file(gemmaSNPInputPath, open = "r")
	
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
	prefix = NULL,
	filepathsPhenoTablePath = NULL,
	annotationPath = NULL){
			
	logRegOutputPath<-snpAnalysisLogReg(
		gemmaGenInputPath = gemmaGenInputPath,
		gemmaSNPInputPath = gemmaSNPInputPath,
		prefix = prefix,
		filepathsPhenoTablePath = filepathsPhenoTablePath
	)



	# Read logistic regression output
	logreg=read.table(logRegOutputPath,header=T, as.is=T)

	position = logreg$position
	pvalue = as.numeric(logreg$negLog10PValue)
	
	manhattan.plot(
		position = position, 
		pvalue = pvalue,
		file = paste(prefix,"_SNP_manhattan.png",sep=""),
		plot.title = paste(prefix," SNP GWAS logistic regression manhatten plot",sep="")
	)	
	
	annot=read.table(file = annotationPath, header=T, sep="\t", quote="")	
	logregSNPAnalysisInfo<-cbind(SNPid = logreg$id,
		position = position,
		negLog10PValue = pvalue,
		ctrl0 = logreg$ctrl0,
		ctrl1 = logreg$ctrl1,
		case0 = logreg$case0,
		case1 = logreg$case1,
		odds = logreg$oddsRatio,
		annot
	)
	
	pValOrder = order(logregSNPAnalysisInfo$negLog10PValue,decreasing=T)
	write.table(logregSNPAnalysisInfo[pValOrder,], file=paste(prefix,"_logreg_annot.txt",sep=""), quote=F, row=F, sep="\t")
	
}

args = commandArgs(trailingOnly = TRUE)
message(length(args))
if(!(length(args) == 6 | length(args) ==0)) {
	stop("\nIncorrect usage\n")
}

if(length(args) == 6){
	source(args[6])
	runLogReg(
		gemmaGenInputPath = args[1],
		gemmaSNPInputPath = args[2],
		prefix = args[3],
		filepathsPhenoTablePath = args[4],
		annotationPath = args[5]
	)

}

