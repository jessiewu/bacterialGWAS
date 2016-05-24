help=paste("",
           "createGEMMAInput.R",
           "Usage: Rscript createGEMMAInput.R [cfPosPath] [bipinfoPath] [imputed.bipinfo] [imputed.bip.patterns] [prefix]",
           "       E.g. Rscript createGEMMAInput.R bact.position_cross_reference.txt bact.bipinfo.txt bact.imputed.bipinfo.txt bact.imputed.bip.patterns.txt bact",
           "",
           "Options:",
           "[cfPosPath]",
           " An output from ClonalFrameML which contains a vector of comma delimited integers.",
           "[bipinfoPath]",
           " An output from snpCalling which contains columns with the headers Position, Allele0, Allele1, A, C, G and T.",
           "[imputed.bipinfo]",
           " An output from compact_SNPs.R which contains columns with the headers Position, Allele0, Allele1, A, C, G and T and Pattern.",
           "[imputed.bip.patterns]",
           " An output from compact_SNPs.R which contains a binary matrix.",
           "[prefix]",
           " File name prefix for output files",
           sep="\n")


createGEMMAInput = function(cfPosPath = NULL, bipinfoPath = NULL, imputed.bipinfo = NULL, imputed.bip.patterns = NULL, prefix = NULL){
	gen.output.file <- paste0(prefix, "_gemma_kinship_genfile.txt")
	snp.output.file <- paste0(prefix, "_gemma_kinship_snpfile.txt")
		
	cf.pos = scan(file=cfPosPath, what=numeric(0), sep=",")
	message("Read in CFML cross reference.")

	bipinfo.df = read.table(file=bipinfoPath, header=T, as.is=T)
	message("Read in bipinfo file before imputation.")
	
	imputed.bipinfo.df = read.table(file=imputed.bipinfo, header=T, as.is=T)
	message("Read in bipinfo file after imputation.")
	
	imputed.bip.patterns.df = read.table(imputed.bip.patterns, header=T, as.is=T)
	message("Read in bip patterns file before imputation.")
	
	cfmlBipPatIndex = cf.pos[bipinfo.df$Position]
	imputedIndex = match(cfmlBipPatIndex,imputed.bipinfo.df$Position)
	bippatIndex = imputed.bipinfo.df$Pattern[imputedIndex]

	gen.file <- imputed.bip.patterns.df[bippatIndex,]
	gen.file <- cbind(paste0("rs",1:nrow(gen.file)), rep(1,nrow(gen.file)), rep(0,nrow(gen.file)), gen.file)
	snp.file <- cbind(paste0("rs",1:nrow(gen.file)), bipinfo.df$Position, rep(24,nrow(gen.file)))
	
	message("Generating GEMMA files.")
	write.table(gen.file, file = gen.output.file, row=F, col=F, sep="\t", quote=F)
	message("Generated GEMMA gen file.")
	write.table(snp.file, file = snp.output.file, row=F, col=F, sep="\t", quote=F)
	message("Generated GEMMA snp file.")
	
		
}

#######################################################################################
## Read in command line arguments
#######################################################################################

args = commandArgs(trailingOnly = TRUE)
if(length(args != 0)){
  if(args[1]=="-help" | args[1]=="-h"){
	  cat(help, sep="\n")
	  q("no")
  }
}

if(length(args) !=5 | length(args)==0) {
	cat(help,sep="\n")
	stop("\nIncorrect usage\n")
}

createGEMMAInput(
	cfPosPath = args[1], 
	bipinfoPath = args[2], 
	imputed.bipinfo = args[3], 
	imputed.bip.patterns = args[4],
	prefix = args[5])