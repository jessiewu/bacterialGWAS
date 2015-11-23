message("run tri/tetra")
help=paste("","tri_tetra_allelic_logreg.R",
"Usage: Rscript tri_tetra_allelic_logreg.R cfml_pos_cross_ref_path cfml_seq_fasta_path prefix filePathPhenoTabPath snpinfofile snppatternsfile",
sep="\n")


### 1. Read the command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args)!=7) {
	cat(help,sep="\n")
	stop("\nIncorrect usage\n")
}



cfml_pos_cross_ref_path = args[1]
cfml_seq_fasta_path = args[2]
prefix=args[3]
filePathPhenoTabPath=args[4]
snpinfo=args[5]
snppatterns=args[6]
annotationFilePath=args[7]




cat(cfml_pos_cross_ref_path,sep="\n")
cat(cfml_seq_fasta_path,sep="\n")
cat(prefix,sep="\n")
cat(filePathPhenoTabPath,sep="\n")
cat(snpinfo,sep="\n")
cat(snppatterns,sep="\n")
cat(annotationFilePath,sep="\n")

# Read txt file
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


##############################################################
## This function is modified for tri/tetra-allelic sites (CHW) 
##############################################################
# Do logistic regression
logreg=function(data,pheno){
	## Fit the models
	fit1 = glm((pheno==1)~1,family=binomial(logit))
	fit2 = glm((pheno==1)~factor(data),family=binomial(logit))

	## Compute the deviance
	dev = abs(fit1$deviance - fit2$deviance)
	## Compute the df
	df = abs(fit1$df.residual - fit2$df.residual)
	
	## Get P-value
	logchi = -log10(pchisq(dev,df,low=F))
	return(logchi)
}



filePathPhenoTab<-read.table(file=filePathPhenoTabPath, header=T)
#filepaths=scan(paste0(filepaths),what=character(0),sep="\n")
filepaths=filePathPhenoTab$filePath
pheno=as.numeric(filePathPhenoTab$phenotype)

##Collects biallelic data?
snpinfo=read.table(snpinfo,h=T,sep="\t",as.is=T)
SNPposition=snpinfo[,1]

##Read in bi-allelic patterns
data = read.table(snppatterns,h=T,sep="\t")

#pheno=as.numeric(scan(pheno,what=character(0),sep="\n"))


cfml_cref=scan(cfml_pos_cross_ref_path,what=integer(0),sep=",")
#cfml=unlist(strsplit(cfml_cref,","))

cfml_imp=read_txt(cfml_seq_fasta_path)
cfml_imp=cfml_imp[1:length(filepaths),]

# Do it in a loop
logreg_output_filename=paste0(prefix,"_tritetra_logreg_output.txt")
logreg_output_log_filename=paste0(prefix,"_tritetra_logreg_output_log.txt")


outputFileHeader = c("id","position",
"ctrl0","ctrl1","ctrl2","ctrl3",
"case0","case1","case2","case3",
"snpPair1","snpPair2","snpPair3","snpPair4","snpPair5","snpPair6",
"oddsRatio1","oddsRatio2","oddsRatio3","oddsRatio4","oddsRatio5","oddsRatio6",
"negLog10PValue")
cat(paste(outputFileHeader, collapse="\t"), 
	file = logreg_output_filename, 
	append=F, 
	sep="\n"
)

start.time<-proc.time()
for(i in 1:nrow(data)) {
	#snpdata = rep(0,ncol(data))
	position=SNPposition[i]
	
	#Read in the allelic variants
	snpinfo2=unlist(snpinfo[i,2:5])
	
	pos_ref=cfml[position]
	data_i = data[i,]
	cfml_i = cfml_imp[,pos_ref] #Get the pattern from CFML output
	
	#Identify the indices with missing values
	whichN = which(data_i==-1)
	if(length(whichN)!=0){
		
		cfml_toimp=cfml_i[whichN]
		
		#############################################
        ## Modified for tri/tetra-allelic sites (CHW)
        #############################################
        for(j in 1:length(snpinfo2)){
        	##Allele index starts from 0.
			cfml_toimp[which(cfml_toimp == snpinfo2[j])] = j - 1 		
		}			
		
		data_i[whichN] = cfml_toimp
	
	}
		
	## Create 2 x J contingency table with pheno on the vertical axis
	tb = table(pheno,as.numeric(data_i))
	control = matrix(0, ncol=4,nrow=1)
	colnames(control) = c(0,1,2,3)
	control[,1:ncol(tb)] = tb[1,]
	cases = matrix(0,ncol=4,nrow=1)
	colnames(cases) = c(0,1,2,3)
	cases[,1:ncol(tb)] = tb[2,]


	#############################################
    ## Modified for tri/tetra-allelic sites (CHW) 
    #############################################
    
	
	## Get alleles from the contingency table,which may not be ordered by allele names.
	## Assuming that the alleles are encoded as 0, 1, 2 ... etc, and they correspond to the
	## columns Allele0, Allele1, Allele2 ... generated from ClonalFrameML. These columns
	## are extracted and stored in the vector snpinfo2.
	alleles<-as.numeric(colnames(tb))
	alleles<-snpinfo2[alleles+1]

	k = 1
	J = ncol(tb)
		
	## Records the pairs for which the odds ratios are calculated.
	#allelePairs = vector(length=choose(J,2))
	allelePairs = rep(NA,6)
	
	## Records the odds-ratios. 
	## The kth odds-ratio in the vector oddsRatios is for 
	## the kth allele-pair in the vector allelePairs.
	#oddsRatios = vector(length=choose(J,2))
		
	oddsRatios = rep(NA,6)	
	Jm1 = J - 1
	
	for(iAllele1 in 1:Jm1){
			
		iAllele2Start = iAllele1 + 1
		for(iAllele2 in iAllele2Start:J){
				
			oddsRatios[k] = (tb[2, iAllele2]*tb[1, iAllele1])/(tb[1, iAllele2]*tb[2, iAllele1])
			##If there are 0s in both rows of the 2 x 2 sub-table
			if(is.nan(oddsRatios[k])){
				oddsRatios[k] = NA
			}
			
			allelePairs[k] = paste0(alleles[iAllele1], alleles[iAllele2])
			k = k + 1
		}
	}

	 
	
	
	info3=c(paste("SNP",i,sep=""), position)
	
	logregression=logreg(unlist(data_i), pheno)
	
	
	#############################################
    ## Modified for tri/tetra-allelic sites (CHW)
    #############################################
    cat(paste( c(info3,control,cases,allelePairs,oddsRatios,logregression), collapse="\t"), 
    	file = logreg_output_filename, append=T, sep="\n")
	
	if((i%%1000)==0) cat("Done",i,"\n")
}
end.time<-proc.time()
diff.time<-end.time-start.time

cat(paste( names(diff.time), collapse="\t"), 
	file = logreg_output_log_filename, append=FALSE, sep="\n")
cat(paste( diff.time, collapse="\t"), 
	file = logreg_output_log_filename, append=TRUE, sep="\n")

annot=read.table(file=annotationFilePath, header=T, sep="\t", quote="")	
logreg.df = read.table(file=logreg_output_filename, header=T, as.is = F)


pval=logreg.df$negLog10PValue
logreg.df= cbind(logreg.df,annot)

o = order(as.numeric(pval), decreasing=T)
write.table(logreg.df[o,], file=paste0(prefix,"_tritetra_logreg_annot.txt"), quote=F, row=F, sep="\t")

position=logreg.df$position

## Plot manhattan plot for log reg
png(paste0(prefix,"_tritetra_SNP_manhattan.png"),width=1000,height=700)
plot(x = position,
	 y = pval, 
	 type="p",
	 main=paste0(prefix," tri and tetra allelic SNP GWAS logistic regression manhatten plot"),
	 xlab="Position in reference genome",
	 ylab="-log10(p-value)"
)
dev.off()