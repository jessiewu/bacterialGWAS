### 1. Read the command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args)!=7) {
		stop("\nIncorrect usage\n")
}



cfmlPath = args[1]
phylogenyPath = args[2]
cfmlFastaFilePath = args[3]
kappa = args[4]
dataFile = args[5]
getMissingSiteCountPath = args[6]
prefix = args[7]

data.df = read.table(file = dataFile, header = T, as.is = T)
fastaFilesTxtPath = paste(prefix,"FastaFilePaths.txt", sep="")
write.table(data.df$filePath, file = fastaFilesTxtPath, row.names=F, col.names = F, quote = F)

noCallSitesPath = paste(prefix, "NoCallSites.txt", sep="")
missingSitesPath = paste(prefix, "MissingSites.txt", sep="")
message("Summary missing site counts started.")
getMissingSiteCountCommand = paste(c("Rscript", getMissingSiteCountPath, fastaFilesTxtPath, noCallSitesPath, missingSitesPath), collapse=" ")
system(getMissingSiteCountCommand)

##TODO: change so it doesn't assume the file name.
totalNoCallSiteCount = scan(file = "totalNoCallSiteCount.txt", what=integer(0), sep="\n")

if(totalNoCallSiteCount == 0){
  cfmlCommand = paste(
	  c(cfmlPath, 
	  phylogenyPath, 
	  cfmlFastaFilePath, 
	  prefix, 
	  "-imputation_only true", 
	  "-kappa",kappa,
	  "-use_incompatible_sites true"),
	  collapse = " "
  )
}else{
  cfmlCommand = paste(
    c(cfmlPath, 
      phylogenyPath, 
      cfmlFastaFilePath, 
      prefix, 
      "-imputation_only true", 
      "-kappa",kappa,
      "-ignore_user_sites", noCallSitesPath,
      "-use_incompatible_sites true"),
    collapse = " "
  )
  
}

system(cfmlCommand)
message("Imputation completed with ClonalFrameML.")
