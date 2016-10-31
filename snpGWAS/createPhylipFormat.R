read_reference_gz = function(ref_file) {
	r = scan(gzfile(ref_file),what=character(0),sep="\n")
	rcat = paste(r[2:length(r)],collapse="")
	return(unlist(strsplit(rcat,"")))
}
read_reference = function(ref_file) {
	r = scan(ref_file,what=character(0),sep="\n")
	rcat = paste(r[2:length(r)],collapse="")
	return(unlist(strsplit(rcat,"")))
}

### 1. Read the command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args)!=5) {
	cat(help,sep="\n")
    stop("\nIncorrect usage\n")
}

reference = args[1]
dataFile = args[2]
bipInfoFile = args[3]
snpInfoFile = args[4]
prefix = args[5]

reference=read_reference(reference)
data.df = read.table(file = dataFile, header = T, as.is = T)
filepaths = data.df$filePath
short_comid = data.df$id

# Variable sites
bipInfo.df = read.table(file=bipInfoFile, h=T, sep="\t", as.is=T)
snpInfo.df = read.table(file=snpInfoFile, h=T, sep="\t", as.is=T)


biallelicSites = bipInfo.df$Position
triTetrallelicSites = snpInfo.df$Position

# Non variable sites
invariantSites = setdiff(1:length(reference), c(biallelicSites, triTetrallelicSites))

# Construct a FASTA and PHYLIP file on the fly
phylip_outfile = paste0(prefix, ".phylip", sep="")

# Read the mapcall file
fa = toupper(read_reference_gz(filepaths[1]))

# Revert non-variable positions to the reference state
fa[invariantSites] = reference[invariantSites]
fa = paste(fa, collapse="")

# Write to PHYLIP
cat(paste(length(filepaths),length(reference)), sep="\n", file = phylip_outfile)
cat(paste(short_comid[1]," ",fa,sep=""), sep="\n", file = phylip_outfile, append=TRUE)
for(i in 2:length(filepaths)) {
	# Read the mapcall file
	fa = toupper(read_reference_gz(filepaths[i]))
	
	# Revert non-core positions to the reference state
	fa[invariantSites] = reference[invariantSites]
	fa = paste(fa, collapse="")
	
	# Write to PHYLIP
	cat(paste(short_comid[i]," ",fa,sep=""),sep="\n",file = phylip_outfile, append=TRUE)	
	cat("Done",i,"\n")
}