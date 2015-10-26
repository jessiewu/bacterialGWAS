# Assumes a fasta file with a single sequence
read_reference = function(ref_file) {
	r = scan(ref_file,what=character(0),sep="\n")
	rcat = paste(r[2:length(r)],collapse="")
	return(unlist(strsplit(rcat,"")))
}

read_reference_gz = function(ref_file) {
	message(ref_file)
	unzipped<-gzfile(ref_file)
	r = scan(unzipped, what=character(0), sep="\n")
	rcat = paste(r[2:length(r)], collapse="")
	lst<-unlist(strsplit(rcat, ""))
	close(unzipped)
	return(lst)
}
### 1. Read the command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args)!=3) {
		    cat(help,sep="\n")
		    stop("\nIncorrect usage\n")
}

prefix = args[1]
reference = args[2]
dataFile = args[3]


reference=read_reference(reference)
data.df = read.table(file = dataFile, header = T, as.is = T)
filepaths=data.df$filePath
short_comid=data.df$id

# First pass count the number of A,C,G,T per site
m = matrix(0,6,length(reference))
for(i in 1:length(filepaths)) {
	# Read the mapcall file
	fa = toupper(read_reference_gz(filepaths[i]))
	m[1,fa=="A"] = m[1,fa=="A"]+1	
	m[2,fa=="C"] = m[2,fa=="C"]+1	
	m[3,fa=="G"] = m[3,fa=="G"]+1	
	m[4,fa=="T"] = m[4,fa=="T"]+1	
	m[5,fa=="N"] = m[5,fa=="N"]+1
	m[6,fa=="-"] = m[6,fa=="-"]+1
	cat("Done",i,"\n")
	
	}

nalleles = (m[1,]>0) + (m[2,]>0) + (m[3,]>0) + (m[4,]>0)
isbiallelic = nalleles==2
is.poly = nalleles>1
is.fixed = !is.poly

# Order the alleles at SNPs
allele.id = matrix(c("A","C","G","T")[apply(m[1:4,is.poly],2,order,decreasing=TRUE)],nrow=4)

# Output filenames
# Files referring to positions in reference:
bip_outfile = paste0(prefix,".bip.patterns.txt");		# biallelic polymorphisms encoded -1 (missing) 0 (allele 0) 1 (allele 1)
bipinfo_outfile = paste0(prefix,".bipinfo.txt");		# positional and allelic information for biallelic polymorphisms
snp_outfile = paste0(prefix,".snp.patterns.txt");		# non-biallelic polymorphisms encoded similarly
snpinfo_outfile = paste0(prefix,".snpinfo.txt");		# non-biallelic polymorphism positional information


# Allocate memory for bip and snp, because need to transform so cannot output on the fly
bip = matrix(NA,sum(isbiallelic),length(short_comid))
snp = matrix(NA,sum(is.poly & !isbiallelic),length(short_comid))
colnames(bip) = short_comid
colnames(snp) = short_comid

# Second pass over the files: populate bip and snp objects
for(i in 1:length(filepaths)) {
	message(filepaths[i])
	# Read the mapcall file
	fa = toupper(read_reference_gz(filepaths[i]))
	fa.poly = fa[is.poly]
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
bip.pat = factor(apply(bip,1,paste,collapse=" "))
snp.pat = factor(apply(snp,1,paste,collapse=""))

# Record only unique patterns, and record the pattern equivalence in the bipinfo file
bip.pat1 = match(levels(bip.pat),bip.pat)
snp.pat1 = match(levels(snp.pat),snp.pat)

# Output compacted bip and snp objects
write.table(bip,bip_outfile,row=FALSE,col=TRUE,quote=FALSE,sep="\t")
write.table(snp,snp_outfile,row=FALSE,col=TRUE,quote=FALSE,sep="\t")

# Output info files
bipinfo = cbind("Position"=which(isbiallelic),"Allele0"=allele.id[1,isbiallelic[is.poly]],"Allele1"=allele.id[2,isbiallelic[is.poly]],"A"=m[1,isbiallelic],"C"=m[2,isbiallelic],"G"=m[3,isbiallelic],"T"=m[4,isbiallelic],"Pattern"=as.numeric(bip.pat)); 
write.table(bipinfo,bipinfo_outfile,row=FALSE,col=TRUE,quote=FALSE,sep="\t")
snpinfo = cbind("Position"=which(is.poly & !isbiallelic),"Allele0"=allele.id[1,!isbiallelic[is.poly]],"Allele1"=allele.id[2,!isbiallelic[is.poly]],"Allele2"=allele.id[3,!isbiallelic[is.poly]],"Allele3"=allele.id[4,!isbiallelic[is.poly]],"A"=m[1,is.poly & !isbiallelic],"C"=m[2,is.poly & !isbiallelic],"G"=m[3,is.poly & !isbiallelic],"T"=m[4,is.poly & !isbiallelic],"Pattern"=as.numeric(snp.pat)); 
write.table(snpinfo,snpinfo_outfile,row=FALSE,col=TRUE,quote=FALSE,sep="\t")

