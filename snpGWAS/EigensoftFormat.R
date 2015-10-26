help=paste("","Eigensoft_format.R takes a file of biallelic SNPs and their positions and formats the data to perform a PCA analysis in Eigensoft","",sep="\n")


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

args = commandArgs(trailingOnly = TRUE)
if(length(args!=0)){
if(args[1]=="-help" | args[1]=="-h"){
	cat(help,sep="\n")
	q("no")
}
}

if(length(args)!=6 | length(args)==0) {
		    cat(help,sep="\n")
		    stop("\nIncorrect usage\n")
}

eigenSoftPath=args[1]
data=args[2]
position=args[3]
prefix=args[4]
fastaPhenoPath = args[5]
npc=args[6]

SNPposition = read.table(position,as.is=T,header=F)[,2]
data = read.table(data,header=F,as.is=T); data=data[,4:ncol(data)]
fastaPhenoPath.df = read.table(file = fastaPhenoPath, header=T,as.is=T)
phenotype=fastaPhenoPath.df$phenotype
fastafiles=fastaPhenoPath.df$filePath



gen_output_filename = paste0(prefix,"_gen.eigenstratgeno")

snp_output_filename = paste0(prefix,"_eigensoftsnp.snp")

ind_output_filename = paste0(prefix,"_eigensoftindiv.ind")


for(i in 1:nrow(data)) {
	datai=as.numeric(data[i,])
	snpdata = rep(0,ncol(data))
	ftr = factor(datai)
	lev = levels(ftr)
	snpdata_allele1 = lev[1]
	snpdata_allele2 = lev[2]
	snpdata[ftr==snpdata_allele1] = 1
	#snpdata[ftr==snpdata_allele2] = 0
	info = c(paste("SNP",i,sep=""), c(24), c(0.0), SNPposition[i])
	
cat(paste(snpdata,collapse=""),file=gen_output_filename,append=(i>1),sep="\n")
	
cat(paste(c(info, snpdata_allele1, snpdata_allele2),collapse="\t"),file=snp_output_filename,append=(i>1),sep="\n")
	
	if((i%%1000)==0) cat("Done",i,"\n")

}


#phenotype=read.table(phenotype,header=F,as.is=T)
phenotype[phenotype==1]="Case"
phenotype[phenotype==0]="Control"

fasta=fastafiles
short_comid=fastaPhenoPath.df$id

ind=cbind(short_comid,c("M"),phenotype)


write.table(ind,file=ind_output_filename,sep="\t",row=F,col=F,quote=F)
message("flag2")
# Construct par file
par=paste(paste0("genotypename: ",getwd(),"/",gen_output_filename),paste0("snpname: ",getwd(),"/",snp_output_filename),paste0("indivname: ",getwd(),"/",ind_output_filename),paste0("evecoutname: ",getwd(),"/",prefix,".pca.evec"),paste0("evaloutname: ",getwd(),"/",prefix,".pca.eval"),paste0("numoutevec: ",npc),paste0("usenorm: NO"),paste0("outliermode: 2"),paste0("numoutlieriter: 0"),paste0("numchrom: 24"),paste0("chrom: 24"),sep="\n")


parFileName = paste(prefix,".parfile.txt", sep="")
cat(par, file = parFileName,sep="\n")
pcaLogFileName = paste(prefix,".pca.logfile", sep="")

system(paste(c(eigenSoftPath,"-p", parFileName, ">", pcaLogFileName), collapse=" "))


