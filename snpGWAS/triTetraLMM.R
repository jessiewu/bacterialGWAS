### 1. Read the command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args)<5 | args[1]=="-h" | args[1]=="-help") {
	cat(help,sep="\n")
	if(args[1]=="-h" | args[1]=="-help"){
	    stop()
	} else {
		 stop("\nIncorrect usage\n")
	}
}

help=paste("triTetraLMM.R performs LMM using GEMMA on tri and tetra allelic sites using a kinship matrix determined by biallelic sites","",
"triTetraPatterns: file ending .snp.patterns produced by the R script compact_SNPs.R - unique SNP patterns",
"kinshipAll: kinship matrix determined by biallelic SNPs",
"phenotypeFiles: file consisting of 7 columns, with the number of rows being the number of phenotypes to test. Column 1: path to phenotype file. Column 2: output file prefix for that phenotype. Column 3: ML log likelihood of the null found in the GEMMA log file when run on the biallelic SNPs. Column 4: snpinfo file from running compact_SNPs. Column 5:snpinfo file on all SNPs. Column 6: cross reference file from running ClonalFrameML. Column 7: snp annotation file.",
"prefix: prefix for gemma formatted files used to test all phenotypes",
"gemmaPath: path to the Wilson modified GEMMA",
"",
"Optional (useful if wanting to run in parallel, default is to test all SNPs):","startLMM: which pattern to begin testing on","endLMM: last pattern to test",
sep="\n")

triTetraPatterns=args[1]
kinshipAll=args[2]
phenotypeFiles=args[3]
#L0=as.numeric(args[4])
prefix=args[4]
gemmaPath=args[5]
if(length(args)>5){
	startLMM=as.numeric(args[6])
	endLMM=as.numeric(args[7])
}



cat(triTetraPatterns,"\n")
cat(kinshipAll,"\n")
cat(phenotypeFiles,"\n")
#cat(L0,"\n")
cat(prefix,"\n")


####################################################################################
## Functions
####################################################################################

getLambda=function(datafiles){
	a=scan(datafiles,what=character(0),sep="\n")
	b=unlist(strsplit(a[13]," "))
	b=b[length(b)]
	return(b)
}

getLogLik=function(datafiles){
	a=scan(datafiles,what=character(0),sep="\n")
	#b=unlist(strsplit(a[13]," "))
	#b=b[length(b)]
	c=unlist(strsplit(a[17]," "))
	c=c[length(c)]
	return(c)
}

getpval=function(datafiles){
	a=read.table(datafiles,header=T)
	a=a$p_lrt
	a=qchisq(a,1,low=F)
	return(a)
}

getLH1=function(datafiles){
	a=read.table(datafiles,header=T)
	a=a$logl_H1
	return(a)
}

####################################################################################
## Check files exist
####################################################################################


file.names=c("triTetraPatterns","kinshipAll","phenotypeFile")
files=c(triTetraPatterns, kinshipAll, phenotypeFiles)

for(i in 1:length(files)){
	if(file.exists(files[i])==FALSE){
		stop(paste0("\nIncorrect usage: ",file.names[i]," doesn't exist\n"))
	}
}


####################################################################################
## Read in files
####################################################################################

triTetraPatterns=read.table(triTetraPatterns,header=T)
#y=scan(phenotypeFile)
phenotypeFiles=read.table(phenotypeFiles,header=F)
pheno_prefix=as.character(phenotypeFiles[,2])
snpinfo=as.character(phenotypeFiles[,4])
snpinfo.full=as.character(phenotypeFiles[,5])
cref=as.character(phenotypeFiles[,6])
annot=as.character(phenotypeFiles[,7])
L0=as.numeric(phenotypeFiles[,3])
phenotypeFiles=as.character(phenotypeFiles[,1])
cat(phenotypeFiles, sep="\n")

for(i in 1:length(phenotypeFiles)){
	if(file.exists(phenotypeFiles[i])==FALSE){
		stop(paste0("\nIncorrect usage: ",pheno_prefix[i]," phenotype file doesn't exist\n"))
	}
}


# Convert patterns to being numeric
triTetraPatterns =apply(triTetraPatterns,1,function(data)as.numeric(data))
# Above function transposes matrix, so need to get back to original dimensions
triTetraPatterns=t(triTetraPatterns)

numAlleles=apply(triTetraPatterns,1,function(data)length(unique(data)))

convertTriTetra=function(data){
#data=as.numeric(data)
if(length(unique(data))==3){
	dataCOV=data
	dataCOV[which(data==1)]=0
	dataCOV[which(data==2)]=1
	data[which(data==2)]=0
} else if(length(unique(data))==4){
	dataCOV=matrix(0,ncol=2,nrow=length(data))
	dataCOV[which(data==2),1]=1
	dataCOV[which(data==3),2]=1
	data[which(data==2 | data==3)]=0
}
return(cbind(data,dataCOV))
}

triTetraPatterns =apply(triTetraPatterns,1,convertTriTetra)

if(length(args) < 6){
 	startLMM=1
 	endLMM=length(triTetraPatterns)
   
}

cat(startLMM,"\n")
cat(endLMM,"\n")

#######################################################################################
## GEMMA input filenames
#######################################################################################


gen_output_filename=paste0(prefix,"_gen_format")
snp_output_filename=paste0(prefix,"_snp_format")
cov_output_filename=paste0(prefix,"_cov_format")


#for(i in 1:length(triTetraPatterns)){
system.time((
for(i in startLMM:endLMM){

					
	# Get info to add for GEMMA gen file
	genInfo =cbind(paste0("triTetra",i),1,0)
				
		
	# create rs numbers for each individual kmerpattern tested
	# make all unique to reduce confusion when combining each of the 'chromosomes'
				
	triTetraInfo =cbind(paste0("triTetra",i),i,24)
				
		
	triTetraPatterns_snpi=matrix(c(genInfo, triTetraPatterns[[i]][,1]),nrow=1)
				
	write.table(triTetraPatterns_snpi,file=gen_output_filename,sep="\t",row=F,col=F,quote=F)
	write.table(triTetraInfo,file=snp_output_filename,sep="\t",row=F,col=F,quote=F)
	write.table(cbind(rep(1,length(triTetraPatterns[[1]][,1])),triTetraPatterns[[i]][,2:ncol(triTetraPatterns[[i]])]),file=cov_output_filename,sep="\t",row=F,col=F,quote=F)
	
	for(j in 1:length(phenotypeFiles)){
		message("dataFile: ", phenotypeFiles[j])
		data.df = read.table(file=phenotypeFiles[j], header=T, as.is=T)
		phenotypeFilePath = paste(c(prefix, "_",pheno_prefix[j],"_phenotype.txt"), collapse="")
		write(data.df$phenotype, file=phenotypeFilePath,  ncolumns=1)
		message("phenotypeFilePath: ", phenotypeFilePath)
		system(paste0(gemmaPath," -g ",gen_output_filename," -p ", phenotypeFilePath," -a ",snp_output_filename," -c ",cov_output_filename," -k ",kinshipAll," -lmm 4 -o ",prefix,"_",pheno_prefix[j],"_gemma_lmmout_triTetra_SNP",i," -maf 0"))
	}

}

))


for(i in 1:length(phenotypeFiles)){

	logfiles=paste0("output/", prefix, "_", pheno_prefix[i],"_gemma_lmmout_triTetra_SNP",startLMM:endLMM,".log.txt")
	assocfiles=paste0("output/", prefix, "_", pheno_prefix[i],"_gemma_lmmout_triTetra_SNP",startLMM:endLMM,".assoc.txt")


	#logLik=as.numeric(sapply(logfiles,getLogLik,USE.NAMES=FALSE))
	#LRTstat =as.numeric(sapply(assocfiles,getpval,USE.NAMES=FALSE))
	LH1 =as.numeric(sapply(assocfiles,getLH1,USE.NAMES=FALSE))

	pvals=rep(0,length(startLMM:endLMM))

	for(j in 1:length(startLMM:endLMM)){
		D=as.numeric(2*(LH1[j]-L0[i]))
		p=pchisq(D,numAlleles[j]-1,low=F)
		pvals[j]=p
	}
	
	outputTable = cbind(paste0("triTetra",startLMM:endLMM),startLMM:endLMM,pvals)
	outputPath = paste0(prefix,"_",pheno_prefix[i],"_triTetra_pvals_patterns_",startLMM,"_to_",endLMM,".txt")
	message(paste("Output file path: ", outputPath))
	
	message(paste(c("Analysis results for ", pheno_prefix[i], "has been processed"), collapse =" "))
	write.table(outputTable, file=outputPath, row=F, col=F, quote=F, sep="\t")
	message(paste(c("Analysis results has been output to  ", outputPath, "successfully."), collapse =" "))
	
	
	assoc=read.table(outputPath,header=F)
	snpinfo.pheno=read.table(snpinfo[i],header=T,sep="\t")
	snpinfo.full.pheno=read.table(snpinfo.full[i],header=T,sep="\t")
	cref.pheno=scan(cref,what=numeric(0),sep=",")
	annot.pheno=read.table(annot[i],header=T,sep="\t",quote="")
	
	assoc=assoc[snpinfo.pheno$Pattern,]
	cref.snps=cref.pheno[snpinfo.full.pheno$Position]
	matchSNPs=match(cref.snps,snpinfo.pheno$Position)
	assoc=assoc[matchSNPs,]
	negLog10=-log10(as.numeric(assoc[,3]))
	
	assoc=cbind(assoc,negLog10,annot.pheno)
	colnames(assoc)[1:4]=c("pattern_rs","pattern_id","pvals","negLog10")
	
	outputPath = paste0(prefix,"_",pheno_prefix[i],"_triTetra_pvals_allSNPs.txt")
	
	write.table(assoc,file=outputPath,row=F, col=T, quote=F, sep="\t")


}





