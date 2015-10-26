read_reference_gz = function(ref_file) {
  unzipped<-gzfile(ref_file)
  r = scan(unzipped, what=character(0), sep="\n")
  rcat = paste(r[2:length(r)], collapse="")
  lst<-unlist(strsplit(rcat, ""))
  close(unzipped)
  return(lst)
}


read_reference = function(ref_file) {
	r = scan(ref_file,what=character(0),sep="\n")
	rcat = paste(r[2:length(r)],collapse="")
	return(unlist(strsplit(rcat,"")))
}

### 1. Read the command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args) != 5) {
		    cat(help,sep="\n")
		    stop("\nIncorrect usage\n")
}

#raxmlPath = "/dipro/mmm/gorm/v3/ana/Cdif/earle/SNPanalysis/RAxML-7.7.6/raxmlHPC-PTHREADS-SSE3"

reference = args[1]
dataFile = args[2]
prefix = args[3]
raxmlPath = args[4]
useConsensus = args[5]

reference=read_reference(reference)
data.df = read.table(file = dataFile, header = T, as.is=T)
filepaths = data.df$filePath
short_comid = data.df$id


# Variable sites
isbiallelic = read.table(paste0(prefix,".bipinfo.txt"),h=T,sep="\t")
isbiallelic=isbiallelic[,1]
# Non variable sites
nonbiallelic = setdiff(1:length(reference),isbiallelic)

# Construct a FASTA and PHYLIP file on the fly
phylip_outfile = paste0(prefix,"_phylip.phylip")

# Read the mapcall file
message(paste0("First file: ",filepaths[1]))
fa = toupper(read_reference_gz(filepaths[1]))
message("First file read")


# Revert non-variable positions to the reference state
fa[nonbiallelic] = reference[nonbiallelic]
fa = paste(fa,collapse="")
# Write to PHYLIP
cat(paste(length(filepaths),length(reference)),sep="\n",file=phylip_outfile)
cat(paste(short_comid[1]," ",fa,sep=""),sep="\n",file=phylip_outfile,append=TRUE)
for(i in 2:length(filepaths)) {
	# Read the mapcall file
	fa = toupper(read_reference_gz(filepaths[i]))
	# Revert non-core positions to the reference state
	fa[nonbiallelic] = reference[nonbiallelic]
	fa = paste(fa,collapse="")
	# Write to PHYLIP
	cat(paste(short_comid[i]," ",fa,sep=""),sep="\n",file=phylip_outfile,append=TRUE)	
	cat("Done",i,"\n")
}


system(paste0(raxmlPath, " -s ", phylip_outfile, " -n ", prefix, " -m GTRCAT -p 12345 -c 1 -T 2 -D ON -f c -F ON -V", sep=""))


