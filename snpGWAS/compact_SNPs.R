# Based on "SNPcalling.R"
# After constructing an ML tree and having used it to do imputation, compact the
# resultant SNP calls into unique patterns to make subsequent analyses more efficient
# in computation time and memory usage.


args = commandArgs(trailingOnly = TRUE)
if(length(args)!=3 | args[1]=="-h" | args[1]=="-help") {
	cat(help,sep="\n")
	if(args[1]=="-h" | args[1]=="-help"){
	    stop()
	} else {
		 stop("\nIncorrect usage\n")
	}
}

help=paste("compact_SNPs.R",
	"After running Clonal Frame ML, compact the SNPs into unique patterns to make subsequent analyses more efficient",
	"",
	"comid: list of comids in the dataset",
	"cfmlFasta: fasta output file from running Clonal Frame ML",
	"outputPrefix: output prefix for output files",
	sep="\n")

dataFile=args[1]
cfmlFasta=args[2]
outputPrefix=args[3]

data.df = read.table(file=dataFile, header=T, as.is=T)


files=c(cfmlFasta, dataFile)
filenames=c(dataFile,cfmlFasta)
for(i in 1:length(files)){
	if(file.exists(files[i])==FALSE){
		stop(paste0("\nIncorrect usage: ",filenames[i]," file doesn't exist\n"))
	}
}

#
# Do this with respect to a particular phenotype. Ie only load typed genomes.
# Assumes a fasta file with a single sequence
read_reference = function(ref_file) {
	r = scan(ref_file,what=character(0),sep="\n")
	rcat = paste(r[2:length(r)],collapse="")
	return(unlist(strsplit(rcat,"")))
}


# List of comids
comid = data.df$id

# Output file
all.comid = substr(scan(pipe(paste0("sed -n 1~2p ",cfmlFasta)),what=character(0)),2,100)
which.keep = 2*match(comid,all.comid)

# for(i in 1:length(all.comid)){
	# message(all.comid[i])
# }

# for(i in 1:length(comid)){
	# message(comid[i])
# }

# stop()

# First pass count the number of A,C,G,T per site
fa = unlist(strsplit(toupper(scan(pipe(paste("sed '",which.keep[1],"q;d' ",cfmlFasta,sep="",collapse="")),what=character(0))),""))
m = matrix(0,4,length(fa))
for(i in 1:length(which.keep)) {
	# Read the mapcall file
	fa = unlist(strsplit(toupper(scan(pipe(paste("sed '",which.keep[i],"q;d' ",cfmlFasta,sep="",collapse="")),what=character(0))),""))
	m[1,fa=="A"] = m[1,fa=="A"]+1
	m[2,fa=="C"] = m[2,fa=="C"]+1
	m[3,fa=="G"] = m[3,fa=="G"]+1
	m[4,fa=="T"] = m[4,fa=="T"]+1
	cat("Done",i,"of",length(which.keep),"\n")
}

nalleles = (m[1,]>0) + (m[2,]>0) + (m[3,]>0) + (m[4,]>0)
is.poly = nalleles>1
is.fixed = !is.poly
isbiallelic = nalleles==2

# Order the alleles at SNPs
allele.id = matrix(c("A","C","G","T")[apply(m[1:4,is.poly],2,order,decreasing=TRUE)],nrow=4)

# Output filenames

bip_outfile = paste0(outputPrefix,".imputed.bip.patterns.txt");		# biallelic polymorphisms encoded -1 (missing) 0 (allele 0) 1 (allele 1)
bipinfo_outfile = paste0(outputPrefix,".imputed.bipinfo.txt");		# positional and allelic information for biallelic polymorphisms
snp_outfile = paste0(outputPrefix,".imputed.ttp.patterns.txt");	# non-biallelic polymorphisms encoded similarly
snpinfo_outfile = paste0(outputPrefix,".imputed.ttpinfo.txt");		# non-biallelic polymorphism positional information

# Allocate memory for bip and snp, because need to transform so cannot output on the fly
bip = matrix(NA,sum(isbiallelic),length(comid))
snp = matrix(NA,sum(is.poly & !isbiallelic),length(comid))
colnames(bip) = comid
colnames(snp) = comid

# Second pass over the files: populate bip and snp objects
for(i in 1:length(comid)) {
	# Read the mapcall file
	fa = unlist(strsplit(toupper(scan(pipe(paste("sed '",which.keep[i],"q;d' ",cfmlFasta,sep="",collapse="")),what=character(0))),""))
	fa.poly = fa[is.poly]
#    This should not be needed (could insert sanity check):
#    fa[is.poly][fa.poly!="A" & fa.poly!="C" & fa.poly!="G" & fa.poly!="T"] = "N"
	# Populate bip object
	bip[,i] = -1
	bip[fa[isbiallelic]==allele.id[1,isbiallelic[is.poly]],i] = 0
	bip[fa[isbiallelic]==allele.id[2,isbiallelic[is.poly]],i] = 1
	# Populate snp object
	snp[,i] = -1
	snp[fa[is.poly & !isbiallelic]==allele.id[1,!isbiallelic[is.poly]],i] = 0
	snp[fa[is.poly & !isbiallelic]==allele.id[2,!isbiallelic[is.poly]],i] = 1
	snp[fa[is.poly & !isbiallelic]==allele.id[3,!isbiallelic[is.poly]],i] = 2
	snp[fa[is.poly & !isbiallelic]==allele.id[4,!isbiallelic[is.poly]],i] = 3
	cat("Done",i,"\n")
}

# Convert BIP and SNP patterns to factors to identify equivalencies
bip.pat = factor(apply(bip,1,paste,collapse=""))
snp.pat = factor(apply(snp,1,paste,collapse=""))

# Record only unique patterns, and record the pattern equivalence in the bipinfo file
bip.pat1 = match(levels(bip.pat),bip.pat)
snp.pat1 = match(levels(snp.pat),snp.pat)

# Output compacted bip and snp objects
write.table(bip[bip.pat1,],bip_outfile,row=FALSE,col=TRUE,quote=FALSE,sep="\t")
write.table(snp[snp.pat1,],snp_outfile,row=FALSE,col=TRUE,quote=FALSE,sep="\t")

# Output info files
bipinfo = cbind("Position"=which(isbiallelic),"Allele0"=allele.id[1,isbiallelic[is.poly]],"Allele1"=allele.id[2,isbiallelic[is.poly]],"A"=m[1,isbiallelic],"C"=m[2,isbiallelic],"G"=m[3,isbiallelic],"T"=m[4,isbiallelic],"Pattern"=as.numeric(bip.pat)); 
write.table(bipinfo,bipinfo_outfile,row=FALSE,col=TRUE,quote=FALSE,sep="\t")
snpinfo = cbind("Position"=which(is.poly & !isbiallelic),"Allele0"=allele.id[1,!isbiallelic[is.poly]],"Allele1"=allele.id[2,!isbiallelic[is.poly]],"Allele2"=allele.id[3,!isbiallelic[is.poly]],"Allele3"=allele.id[4,!isbiallelic[is.poly]],"A"=m[1,is.poly & !isbiallelic],"C"=m[2,is.poly & !isbiallelic],"G"=m[3,is.poly & !isbiallelic],"T"=m[4,is.poly & !isbiallelic],"Pattern"=as.numeric(snp.pat)); 
write.table(snpinfo,snpinfo_outfile,row=FALSE,col=TRUE,quote=FALSE,sep="\t")


