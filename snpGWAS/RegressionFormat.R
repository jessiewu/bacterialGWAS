###############################################################################
## Read txt file
###############################################################################
read_txt = function(txt_file){
	d = scan(txt_file, what=character(0), sep="\n")
	d = d[seq(2,length(d),by=2)]
	d1 = unlist(strsplit(d[1],""))
	data = matrix("",length(d),length(d1))
	data[1,] = d1
	for(i in 2:length(d)) {
		data[i,] = unlist(strsplit(d[i],""))
		if((i%%100)==0) cat("Done",i,"\n")
	}
	return(data)
}


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
## Do logistic regression
###############################################################################
imputeMissingValues<-function(
	data_i = NULL,
	cfml_i = NULL,
	bipinfo_i = NULL){
		#message(bipinfo_i)
	
	
	#get position of missing values
	whichN = which(data_i==-1)
	
	
	if(length(whichN) > 0){
			
		cfml_toimp = cfml_i[whichN]
		for(j in 1:length(cfml_toimp)){
			if(cfml_toimp[j] == bipinfo_i[1]){
				cfml_toimp[j] = 0
			} else {
				cfml_toimp[j] = 1
			}
		}
	
		data_i[whichN] = cfml_toimp	
	}
	return(data_i)
	
}

###############################################################################
## Do logistic regression
###############################################################################
imputeMultiAllelicMissingValues<-function(
	data_i = NULL,
	cfml_i = NULL,
	snpinfo_i = NULL){
	
	
	#Identify the indices with missing values
	whichN = which(data_i==-1)
	if(length(whichN)!=0){
		
		cfml_toimp = cfml_i[whichN]
		
		#############################################
        ## Modified for tri/tetra-allelic sites (CHW)
        #############################################
        for(j in 1:length(snpinfo_i)){
        	##Allele index starts from 0.
			cfml_toimp[which(cfml_toimp == snpinfo_i[j])] = j - 1 		
		}			
		
		data_i[whichN] = cfml_toimp
	
	}

	return(data_i)	
}





###############################################################################
## Get a vector of binary alleles and identify the major and minor alleles
###############################################################################
getAlleles<-function(data_i = NULL){
	ftr = factor(as.character(data_i))
	lev = levels(ftr)
	snpdata_allele1 = lev[1]
	snpdata_allele2 = lev[2]
	
	countallele1 = length(which(data_i==snpdata_allele1))
	countallele2 = length(which(data_i==snpdata_allele2))
	
	if(countallele1 <= countallele2){
		mallele = snpdata_allele1
		nonmallele = snpdata_allele2
		} else {
			mallele = snpdata_allele2
			nonmallele = snpdata_allele1
		}
	alleles<-list(minor.allele = mallele, major.allele = nonmallele)
	return(alleles)	
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
## Create output for an analysis in GEMMA
###############################################################################
createGemmaOutputFile<-function(
	iteration = NULL,
	snpData = NULL,
	alleles = NULL,
	gen_output_filename = NULL){
	
	info1 = c(paste("rs", iteration, sep=""), alleles$minor.allele, alleles$major.allele)
	
	cat(paste(c(info1, snpData), collapse="\t"), 
		file = gen_output_filename, append = (iteration>1), sep="\n")
}

###############################################################################
## Create output for an analysis in LMM with tri/tetrallelic sites
###############################################################################
createTempTriTetrallelicLMMOutputFile<-function(
	iteration = NULL,
	snpData = NULL,
	tempOutputFilename = NULL){
			
	cat(paste(snpData, collapse="\t"), 
		file = tempOutputFilename, 
		append = (iteration>1), sep="\n"
	)
}

###############################################################################
## Create with all the SNP positions
###############################################################################
createSNPOutputFile<-function(
	iteration = NULL,
	snpPosition = NULL,
	snp_output_filename = NULL){
		
	info2 = c(paste("rs", iteration, sep=""), snpPosition)
	#message(paste("snp_output_filename: ",snp_output_filename))
	cat(paste(c(info2, c(24)),collapse="\t"), 
		file = snp_output_filename, append = (iteration>1), sep="\n")
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
		append = (iteration>1), 
		sep="\n"
	)
	#message("bye")
}

###############################################################################
## Performs logistic regression and outputs summary for a single snp
###############################################################################
createGEMMAInputFormat<-function(
	data_i = NULL,
	cfml_i = NULL,
	bipinfo_i = NULL,
	gen_output_filename = NULL,
	snp_output_filename = NULL,
	logreg_output_filename = NULL,
	snpPosition = NULL,
	iteration = NULL){
		
	data_i = imputeMissingValues(
		data_i = data_i,
		cfml_i = cfml_i,
		bipinfo_i = bipinfo_i
	)
		
	## Create vector
	snpData = rep(0,length(data_i))
	
	alleles<-getAlleles(data_i)	
	ftr = factor(as.character(data_i))
	snpData[ftr==alleles$minor.allele] = 1
	

	createGemmaOutputFile(
		iteration = iteration,
		snpData = snpData,
		alleles = alleles,
		gen_output_filename = gen_output_filename
	)
	
	createSNPOutputFile(
		iteration = iteration, 
		snpPosition = snpPosition, 
		snp_output_filename = snp_output_filename
	)

}


###############################################################################
## Performs logistic regression and outputs summary for a single snp
###############################################################################
createTriTetraAllelicLMMInputFormat<-function(
	data_i = NULL,
	cfml_i = NULL,
	snpinfo_i = NULL,
	snpPosition = NULL,
	tempOutputFilename = NULL,
	iteration = NULL){
		
	snpData = imputeMultiAllelicMissingValues(
		data_i = data_i,
		cfml_i = cfml_i,
		snpinfo_i = snpinfo_i
	)	

	createTempTriTetrallelicLMMOutputFile(
		iteration = iteration,
		snpData = snpData,
		tempOutputFilename  = tempOutputFilename
	)	

}

###############################################################################
## Create the input file in the format required for analyses in GEMMA
###############################################################################
createGEMMAInputFile<-function(
	cfml.pos.cross.ref.path = NULL,
	cfml.ml.seq.fa.path = NULL,
	bipinfo.path = NULL,
	bip.pattern.path = NULL,
	prefix = NULL,
	filepathsPhenoTablePath = NULL){
		
	filepathsPhenoTable<-read.table(file=filepathsPhenoTablePath, header=T)
	#message(names(filepathsPhenoTable))
		
	#Get file paths	
	filepaths=filepathsPhenoTable$filePath
	#Get phenotype
	#phenotype=as.numeric(filepathsPhenoTable$phenotype)
	#Get bipinfo
	bipinfo=read.table(bipinfo.path, header=T, sep="\t")
	SNPposition=bipinfo[,1]
	
	
	data = read.table(bip.pattern.path, header=T, as.is=T)

	cfml=scan(cfml.pos.cross.ref.path, what=integer(0),sep=",")
	cfml_imp=read_txt(cfml.ml.seq.fa.path)
	#message(paste("dim: ",dim(cfml_imp)))
	#message(paste("length: ",length(filepaths)))
	cfml_imp=cfml_imp[1:length(filepaths),]
	#message(paste("dim: ",dim(cfml_imp)))


	gen_output_filename = paste0(prefix,"_gemma_gen_format.txt")
	snp_output_filename = paste0(prefix,"_gemma_snp_format.txt")
	logreg_output_filename = paste0(prefix,"_logreg_output.txt")
	
	snp.count = nrow(bipinfo)
	message(paste("Total number of binary snps: ",snp.count))
	for(i in 1:snp.count) {
		
		## get bip position
		position = SNPposition[i]		
		#get raw data for the ith snp
		data_i = data[i,]
		#get the index position in the cfml imputated patterns
		pos_ref = as.numeric(cfml[SNPposition[i]])
		#message(paste("pos_ref",pos_ref))
		#get replacement for the missing values of the ith snp
		cfml_i = cfml_imp[,pos_ref]	
	
		createGEMMAInputFormat(
			data_i = data_i, 
			#phenotype = phenotype,
			cfml_i = cfml_i, 
			bipinfo_i = bipinfo[i, 2:3],
			gen_output_filename = gen_output_filename,
			snp_output_filename = snp_output_filename,
			logreg_output_filename = logreg_output_filename,
			snpPosition = SNPposition[i],
			iteration = i
		)
	
		if((i%%1000)==0) cat("Done",i,"\n")
	}
	
	outputFilePaths<-list(
		gen_output = gen_output_filename,
		snp_output = snp_output_filename,
		logreg_output = logreg_output_filename
	)
	
	return(outputFilePaths) 
}

###############################################################################
## Create the output file in the format required for 
## LMM analyses with tri/tetra-allelic snps
###############################################################################
createTriTetralleicLMMInputFile<-function(
	cfml.pos.cross.ref.path = NULL,
	cfml.ml.seq.fa.path = NULL,
	snpinfo.path = NULL,
	snp.pattern.path = NULL,
	prefix = NULL,
	filepathsPhenoTablePath = NULL){
		
	filepathsPhenoTable<-read.table(file=filepathsPhenoTablePath, header=T)
		
	#Get file paths	
	filepaths=filepathsPhenoTable$filePath

	#Get snpinfo
	snpinfo = read.table(snpinfo.path, header=T, sep="\t", as.is=T)
	SNPposition = snpinfo$Position
	
	
	data = read.table(snp.pattern.path,h=T,sep="\t")

	cfml=scan(cfml.pos.cross.ref.path, what=integer(0),sep=",")
	cfml_imp=read_txt(cfml.ml.seq.fa.path)
	cfml_imp=cfml_imp[1:length(filepaths),]
	
	lmmFormatTempOutputFilename = paste0(prefix,"_temp_tritetrallelic_LMM_format.txt")
	lmmFormatOutputFilename = paste0(prefix,"_tritetrallelic_LMM_format.txt")
		
	snp.count = nrow(snpinfo)
	message(paste("Total number of tri/tetrallelic snps: ",snp.count))
	allele.headers<-c('Allele0','Allele1','Allele2','Allele3')
	
	for(i in 1:snp.count) {
		
		## get bip position
		position = SNPposition[i]
		#get the index position in the cfml imputated patterns
		pos_ref = as.numeric(cfml[position])		
		#get raw data for the ith snp
		data_i = data[i,]
		#get replacement for the missing values of the ith snp
		cfml_i = cfml_imp[,pos_ref]	
	
		createTriTetraAllelicLMMInputFormat(
			data_i = data_i,
			cfml_i = cfml_i,
			snpinfo_i = unlist(snpinfo[i,allele.headers]),
			snpPosition = position,
			tempOutputFilename = lmmFormatTempOutputFilename,
			iteration = i
		)
	
		if((i%%1000)==0) cat("Done", i, "\n")
	}
	
	
	tempOutput.df<-read.table(file=lmmFormatTempOutputFilename,header=F)
	tempOutput.df<-t(tempOutput.df)
	nrow.temp<-nrow(tempOutput.df)
	nrow.pheno<-nrow(filepathsPhenoTable)
	if(nrow.temp != nrow.pheno){
		message(paste(c("The temporary multi-allelic table and the phenotype table have different lengths:", 
			nrow.temp, nrow.pheno), collapse=" "))
		stop()
	}
	names(tempOutput.df)<-SNPposition
	lmmOutput.df<-cbind(filepathsPhenoTable, tempOutput.df)
	write.table(lmmOutput.df,
		file=lmmFormatOutputFilename,
		sep="\t",
		row.names=FALSE,
		quote=FALSE
	)
	
	
	return(lmmFormatOutputFilename) 
}




###############################################################################
###############################################################################
###############################################################################

### 1. Read the command line arguments
args = commandArgs(trailingOnly = TRUE)
if(!(length(args) ==6 | length(args) ==0)) {
	    stop("\nIncorrect usage\n")
}


start.time = proc.time()
if(length(args) == 6){
	cfml.pos.cross.ref.path = args[1]
	cfml.ml.seq.fa.path = args[2]
	bipinfo.path = args[3]
	bip.pattern.path = args[4]
	prefix = args[5]
	filepathsPhenoTablePath = args[6]

   

	cat(cfml.pos.cross.ref.path,sep="\n")
	cat(cfml.ml.seq.fa.path,sep="\n")
	cat(prefix,sep="\n")
	cat(filepathsPhenoTablePath,sep="\n")


	createGEMMAInputFile(
		cfml.pos.cross.ref.path = cfml.pos.cross.ref.path,
		cfml.ml.seq.fa.path = cfml.ml.seq.fa.path,
		bipinfo.path = bipinfo.path,
		bip.pattern.path = bip.pattern.path,
		prefix = prefix,
		filepathsPhenoTablePath = filepathsPhenoTablePath
	)
}
end.time = proc.time()

diff.time = end.time - start.time

message(diff.time[1])
message(diff.time[2])
message(diff.time[3])




