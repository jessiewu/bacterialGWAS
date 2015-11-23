## SNP logistic regression

help = paste("","do_logreg.R runs a logistic regression and outputs a -log10(pvalue), numbers of cases and controls, and an odds ratio","Usage: Rscript do_logreg.R -bips -phenotype -prefix","       E.g. Rscript do_logreg.R -bips Cdifbips.txt -phenotype pheno.txt ","","Options:","-bips"," File containing biallelic SNPs called as 1 or 0 against a reference. One row per biallelic SNP, 'SNPid, bip bip bip'"," E.g. 'rs1	1	0	0	0	0	0	1	0	0' in a tabular format","-phenotype"," File containing binary phenotype coded as 0 for control and 1 for case","-prefix"," Output file name prefix","",sep="\n")

# Do logistic regression
logreg=function(data,pheno){
		logr1=glm((pheno==1)~1,family=binomial(logit))$deviance
		logr2=glm((pheno==1)~1+data,family=binomial(logit))$deviance
		dev=abs(diff(c(logr1,logr2)))
		logchi = -log10(pchisq(dev,1,low=F))
		return(logchi)
}



### 1. Read the command line arguments
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

s=seq(2,by=2,length(args))
args2=args[seq(1,by=2,length(args))]
inputs=matrix(0,ncol=2,nrow=length(args2))
for(i in 1:length(args2)){
	arg=unlist(strsplit(args2[i],""))
	inputs[i,1]=paste(arg[2:length(arg)],collapse="")
	inputs[i,2]=args[s[i]]
}

bips=c(0)
phenotype=c(0)
prefix=c(0)


for(i in 1:length(inputs[,1])){
	if(inputs[i,1]=="bips"){
		bips=inputs[i,2]
		}
		else if(inputs[i,1]=="phenotype"){
		phenotype=inputs[i,2]
		}
		else if(inputs[i,1]=="prefix"){
		prefix=inputs[i,2]
	}
}

bips_file=bips
existance = file.exists(bips_file)
if(length(which(existance==FALSE))>0) {
		    stop("\nIncorrect usage: Bips file doesn't exist\n")
}
bips=read.table(bips,header=F,sep="\t", as.is=T, quote ="\"")
SNPid=bips[,1]

unique.count<-length(unique(bips[,1]))
bip.count<-length(bips[,1])
bips=bips[,2:ncol(bips)]

existance = file.exists(phenotype)
if(length(which(existance==FALSE))>0) {
		    stop("\nIncorrect usage: phenotype file doesn't exist\n")
}
pheno=as.numeric(scan(phenotype,what=character(0),sep="\n"))

logreg_output_filename=paste0(prefix,"_SNP_logreg_output.txt")

cat(paste0("Bips file: ",bips_file),paste0("Phenotype: ",phenotype),paste0("Prefix: ",prefix),paste0("Start time: ",Sys.time()),file=paste0(prefix,"_input_logfile.txt"),sep="\n")

cat(paste0("Saved input and run information to ",prefix,"_input_logfile.txt"),sep="\n")

headerinfo=c("ID","Controls0","Controls1","Cases0","Cases1","Pvalue","Oddsratio")

cat(paste(c(headerinfo),collapse="\t"),file=logreg_output_filename,append=FALSE,sep="\n")

for(i in 1:nrow(bips)) {
	
	data=as.numeric(bips[i,])
	
	id=SNPid[i]
	
	tb=table(pheno,data)
	controls=tb[1,]
	cases=tb[2,]
	if(length(which(tb[1,]==0)) ==0 | length(which(tb[2,]==0)) ==0){
		odds=(tb[2,2]*tb[1,1])/(tb[1,2]*tb[2,1])
	} else {
		odds="NA"
	}
	
	logregression=logreg(data,pheno)
	
	z<-c(id,controls,cases,logregression,odds)
	if(length(z) != 7){
		stop(z)
	}	
	
	cat(paste(z,collapse="\t"),file=logreg_output_filename,append=TRUE,sep="\n")	
	
	if((i%%1000)==0) cat("Done",i,"\n")
}

cat(paste0("Logistic regression completed: ",Sys.time()),file=paste0(prefix,"_input_logfile.txt"),append=TRUE,sep="\n")


# Read logistic regression output
logreg=read.table(paste0(prefix,"_SNP_logreg_output.txt"),header=T,as.is=TRUE,sep="\t")

SNPid = logreg[,1]
pvalue = as.numeric(logreg[,6])


## Plot manhattan plot for log reg
png(paste0(prefix,"_SNP_manhattan.png"),width=1200,height=750)
plot(x = 1:length(SNPid), y = pvalue, type="p",main=paste0(prefix," GWAS logistic regression manhatten plot"),xlab="SNPid",ylab="-log10(p)")
dev.off()
