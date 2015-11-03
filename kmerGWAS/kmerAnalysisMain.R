# File KmerAnalysisMain.R
# Authors: Earle, S., Wu, C.-H. and Wilson, D. J. 
#
# Copyright (C) 2015 University of Oxford
#
# This file is part of the bacterial GWAS pipeline.
# See the NOTICE file distributed with this work for additional
# information regarding copyright ownership and licensing.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
#  This is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this software package; if not, write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA  02110-1301  USA

help = paste(
  "",
  "KmerAnalysisMain.R performs a kmer GWAS",
  "Authors: Earle, S., Wilson, D. J. and Wu, C.-H.",
  "Usage: Rscript kmerAnalysisMain.R -dataFile -srcDir -prefix -genFileFormat -externalSoftware -refSeqFasta1 -ncbiSummary -db2 -relateMatrix",
  "       E.g. Rscript kmerAnalysisMain.R -dataFile ecol_id_bam_pheno.txt -srcDir /home/scripts/R/kmer/ -prefix ecol -externalSoftware extSoftware.txt -genFileFormat BAM -refSeqFasta1 ecolRefSeq1.txt -ncbiSummary ecolNcbiSummary.txt -db2 /home/ecol/db/ecolDB2 -relateMatrix /home/eco/lmm/output/",
  "","Options:",
  "-dataFile",
  " The path of a file containing list of file paths to bam files or kmer files, corresponding ids and phenotypes.",
  " If the kmer files are to be created, then the paths to the corresponding bam files needs to be provided.",
  " Otherwise the paths to the kmer files needs to be provided.",
  "-srcDir",
  " Path of the directory that contains all the R scripts.",  
  "-genFileFormat",
  "The file format of the genomic data.",
  "-minCov",
  "The minimum depth of kmers (E.g. 1 if based on assemblies).",
  "-prefix",
  " The prefix of the output files.",  
  "-runLMM",
  " Whether or not to run the analysis LMM (accepted value: TRUE or FALSE; default = TRUE)",
  "-externalSoftware",
  " The path of a comma delimited file containing the name and paths of the external software used in the analysis",
  "-createDB",
  " Whether or not to create the databases for annotation (accepted value: TRUE or FALSE; default = TRUE)",
  "-signif",
  " The number of top significant kmers to be annotated (default = 10000)",
  "-kmerFileDir",
  " The path of the directory where the created kmers are to be stored (default = current directory).",
  "-nproc",
  " The number of processors (default = 1).",
  "-refSeqFasta1",
  " The path of the file that lists the paths of the genomes to be used to create the first database for annotation.",
  "-refSeqFasta2",
  " The path of the file that lists the paths of the genomes to be used to create the second database for annotation.",
  "-ncbiSummary",
  " The summary file downloaded from the NCBI database for the genes of a given species",
  "-db1",
  " The path to the first database used for annotation",
  "-db2",
  " The path to the second database used for annotation",
  "-ncbiAnnot",
  " The path to the database created from the NCBI summary file",
  "-relateMatrix",
  " The path of the relateMatrix created from GEMMA",
  "-signifLMM",
  " The number of top significant kmers to be used for the LMM analysis (default = 100000)",
  "",
  sep="\n"
)

################################################################################################
## Format the dir appropriately.
## @dir: The directory path to be formated properly.
################################################################################################
formatDir = function(dir = NULL){
  dirLength = nchar(dir)
  if(substr(dir, dirLength, dirLength) !='/'){
    dir = paste(dir,"/",sep="")
  }
  return(dir)
}


################################################################################################
## Chech the existence of file paths.
## @filePath: The path of the file to be checked.
################################################################################################
checkExistence = function(filePath = NULL){
  message(filePath)
  message(class(filePath))
	doesNotExist = which(!file.exists(filePath))
  message("Done!")
	if(length(doesNotExist) > 0){
		stop(paste("The file", filePath[doesNotExist]," does not exist!", collapse="\n", sep=""))
	}
	
}


##################################################################################
## Turn command line inputs into a n x 2 matrix
## @args: A string vector of command line arguments.
##################################################################################
getCommandLineInputMatrix = function(args = NULL){
  argsCount = length(args)/2
  inputs = matrix(nrow=argsCount, ncol=2)
  inputs[,1] = args[c(1:argsCount)*2 -1]
  
  if(length(which(regexpr("-",inputs[,1]) != 1)) > 0){
    stop("Argument names must start with '-'! E.g. -phylogeny.")
  }
  
  inputs[,1] = gsub("-", "", inputs[,1])
  inputs[,2] = args[c(1:argsCount)*2]
  return(inputs)
  
}


##################################################################################
## Extract input arguments.
## @argName: Name of the argument to be extracted
## @commandLineInputs: A n x 2 matrix, where n is the number of command line inputs.
## @default: default value of the input argument if not provided
## The first column of the matrix contains the names of the arguments, and the second
## column contains the arguments values.
##################################################################################
extractInputArgument = function(argName = NULL, commandLineInputs = NULL,
                                default = NULL, canBeNULL = FALSE){
  
  argIndex = which(commandLineInputs[,1] == argName)
  argIndexCount = length(argIndex)
  if(argIndexCount == 0){
    if(canBeNULL | !is.null(default)){
      return(default)
    }else{
      stop(paste(c("The argument", argName, "must be specified!"), collapse=" "))
    }
  }else if(argIndexCount > 1){
    stop(paste(c("The argument", argName, "has been specified multiple times!"), collapse=" "))
  }else{
    return(commandLineInputs[argIndex,2])
  }
  
}


##################################################################################
## Extract paths of the software to be used.
## @softwareName: The name of the software to be used.
## @pathTB.df: A data frame containing a list of external software and their paths
##################################################################################
getSoftwarePath = function(
  softwareName = NULL, 
  pathTb.df = NULL){	
  
  path = pathTb.df$path[which(toupper(pathTb.df$name) == toupper(softwareName))]
  checkExistence(path)
  
  return(path)
}


##################################################################################
## Extract paths of the software to be used.
## @createDB: The input indicating whether or not to create the databases for annotation
## @refSeqFasta1: Path of a file containing a list of genomes for creating the first database for annotation
## @refSeqFasta2: Path of a file containing a list of genomes for creating the second database for annotation
## @ncbiSummary: Path of the summary file downloaded from NCBI
## @db1: Path prefix of the first database
## @db2: Path prefix of the second database
## @ncbiAnnot: Path of the NCBI database
##################################################################################
checkInputRelatedToDB = function(createDB = NULL,
                                 refSeqFasta1 = NULL,
                                 refSeqFasta2 = NULL,
                                 ncbiSummary = NULL,
                                 db1 = NULL,
                                 db2 = NULL,
                                 ncbiAnnot = NULL){
  
  if(createDB){
    if(is.null(refSeqFasta1) & is.null(db2)){
      stop(paste(c("The list of genomes (refSeqFasta1) for creating the first database must be specified if the databases are to be created,",
                   "or the user needs to provide an existing database as the -db1 argument."), collapse="\n"))
    }else{
      if(is.null(refSeqFasta2) & is.null(db2)){
        
        stop(paste(c("If the databases are to be created (createDB = TRUE), the user needs to provide either:", 
                     "A list of genomes (refSeqFasta2) for creating the database or",
                     "An existing second database (db2)"), collapse="\n"))
        
      }
    }
    
    if(is.null(ncbiSummary)){
      stop("The NCBI summary file must be specified.")
    }
    
  }else{
    dbNames = c("db1", "db2","ncbiAnnot")
    missing = rep(FALSE, 3)
    if(is.null(db1)){
      missing[1] = TRUE
    }
    if(is.null(db2)){
      missing[2] = TRUE
    }
    if(is.null(ncbiAnnot)){
      missing[3] = TRUE
    }
    
    if(any(missing)){
      stop(paste(c("The following arguments must be specified if not creating the databases:", dbNames[missing]), collapse="\n"))
    }  
    
  }
  message("Database related input checks completed.")
}


##########################################################################################
## Create kmer files from bam files.
## @bam_files: A list of bam files of the genomes (one bam file per genome)
## @comid: A list of id that corresponds to the genomes
## @phenotype: A list of phenotypes that corresponds to the list genomes
## @kmerFileDir: The output directory of the kmer files created
## @prefix: The prefix of the output file containing the kmer file paths and phenotype
##########################################################################################
createKmerFilesFromBAM = function(bam_files = NULL, 
                           comid = NULL, 
                           phenotype = NULL,
                           nproc = NULL,
                           kmerFileDir = NULL,
                           prefix = NULL,
                           externalSoftwarePaths = NULL){
          
  list_len = ceiling(length(bam_files)/nproc)
  
  for(i in 1:nproc) {
    bam_proc = bam_files[((i-1)*list_len+1):pmin(length(bam_files),i*list_len)]
    comid_proc = comid[((i-1)*list_len+1):pmin(length(comid),i*list_len)]
    write(bam_proc, paste0(prefix, "_bam_", i, "c.txt"), ncol=1)
    write(comid_proc, paste0(prefix, "_comids_", i, "c.txt"), ncol=1)
  }  
  
  split_bam=paste0(prefix, "_bam_", 1:nproc, "c.txt")
  split_bam=paste(split_bam, collapse = " ")
  split_comids=paste0(prefix, "_comids_", 1:nproc, "c.txt")
  split_comids=paste(split_comids, collapse = " ")
  
  if(nproc > 1){
    parallelInput = paste(c(doAndSortDskR, "-bamFilePaths", split_bam, "-ids", split_comids, "-externalSoftware", externalSoftwarePaths, 
            "-kmerFileDir", kmerFileDir, "-qualityScore", qualityScore), collapse = " ::: ")
    
    system(paste(c(parallelPath, "-j", nproc, "--xapply", "Rscript", parallelInput), collapse = " "))
    
  }else{
    message(externalSoftwarePaths)
    system(paste(c("Rscript", doAndSortDskR, "-bamFilePaths", split_bam[1], "-ids", split_comids[1], "-externalSoftware", externalSoftwarePaths, "-kmerFileDir", kmerFileDir, "-qualityScore", qualityScore), collapse=" "))
    
  }
  
  kmerFilePaths = paste(kmerFileDir, comid, ".kmer31.txt.gz", sep="")
  
  kmerPheno.df = data.frame(filePath=kmerFilePaths, phenotype = phenotype)
  kmerPhenoPath = paste(c(formatDir(getwd()), prefix,"KmerPheno.txt"), collapse="")
  write.table(kmerPheno.df, file = kmerPhenoPath, row.names=F, col.names=T, sep="\t", quote=F)
  
  return(kmerPhenoPath)
  
  
  
}

##########################################################################################
## Create kmer files from bam files.
## @dataFilePath: A list of bam files of the genomes (one bam file per genome)
## @kmerFileDir: The output directory of the kmer files created
## @extSoft: The file that specifies the path of the used external software
##########################################################################################
createKmerFilesFromFASTA = function(createKmersFromFastaR = NULL, dataFilePath = NULL, kmerFileDir = NULL, extSoft = NULL){
	
	system(paste(c("Rscript", createKmersFromFastaR, "-dataFile", dataFilePath, "-externalSoftware", extSoft, "-kmerFileDir", kmerFileDir), collapse=" "))
	datTab.df = read.table(file = dataFilePath, header = T, as.is = T)
	kmerFilePaths = paste(kmerFileDir, datTab.df$id, ".kmer31.txt.gz", sep="")
	
	kmerPheno.df = data.frame(filePath=kmerFilePaths, phenotype = datTab.df$phenotype)
  	kmerPhenoPath = paste(c(formatDir(getwd()), prefix,"KmerPheno.txt"), collapse="")
  	write.table(kmerPheno.df, file = kmerPhenoPath, row.names=F, col.names=T, sep="\t", quote=F)
  
  	return(kmerPhenoPath)

	
}


##########################################################################################
## Run kmer analysis
## @kmerAnalysisR: The path of the R script kmerAnalysis.R
## @kmerPhenoPath: Path of the file containing the kmer files and the binary phenotypes
## @prefix: Prefix of the output files
## @removeTxt: Whether or not to remove the (usually) gigantic kmer file
##########################################################################################
runKmerAnalysis = function(kmerAnalysisR = NULL,
                           kmerPhenoPath = NULL, 
                           prefix = NULL, 
                           removeKmerTxt = NULL,
                           externalSoftwarePaths = NULL,
                           minCov = NULL){
  message("minCov: ", minCov)
  system(paste(c("Rscript", kmerAnalysisR, "-dataFile", kmerPhenoPath, "-prefix", prefix, 
                 "-removeKmerTxt", removeKmerTxt, "-minCov", minCov, "-externalSoftware", externalSoftwarePaths), collapse=" "))
  
  # This is not ideal as this requires knowing the suffix
  # and how the names are generated before hand. It would
  # be better if this is some returned value.
  dir = formatDir(getwd())
  outputList = list(patternKey = paste(dir, prefix, ".gwaskmer-out.patternKey.txt", sep=""),
       patternIndex = paste(dir, prefix, ".gwaskmer-out.patternIndex.txt", sep=""),
       chisqStat = paste(dir, prefix, ".gwaskmer-out.chisqStat.txt", sep=""),
       meanDepthCase = paste(dir, prefix, ".gwaskmer-out.meanDepthCase.txt", sep=""),
       meanDepthCtrl = paste(dir, prefix, ".gwaskmer-out.meanDepthCtrl.txt", sep=""),
       nPresentCase = paste(dir, prefix, ".gwaskmer-out.nPresentCase.txt", sep=""),
       nPresentCtrl = paste(dir, prefix, ".gwaskmer-out.nPresentCtrl.txt", sep=""),
       kmer = paste(dir, prefix, ".gwaskmer-out.kmer.txt", sep=""),
       rgrssLogP = paste(dir, prefix, ".gwaskmer-out.rgrssLogP.txt", sep=""),
       trendLogP = paste(dir, prefix, ".gwaskmer-out.trendLogP.txt", sep=""),
       chisqLogP = paste(dir, prefix, ".gwaskmer-out.chisqLogP.txt", sep=""),
       rgrssStat = paste(dir, prefix, ".gwaskmer-out.rgrssStat.txt", sep=""),
       trendStat = paste(dir, prefix, ".gwaskmer-out.trendStat.txt", sep=""))  
  return(outputList)
  
}


##########################################################################################
## Create databases for kmer annotation
## @parseNcbiGeneR: The path of the R script parse_ncbi_gene.R
## @refSeqFasta1: The path of the fasta files for building the first database
## @refSeqFasta2: The path of the fasta files for building the second database
## @ncbiSummary: The summary file downloaded from ncbi
## @prefix: The prefix of databases to be built
##########################################################################################
runParseNcbiSummary = function(parseNcbiGeneR = NULL,
                               refSeqFasta1 = NULL,
                               refSeqFasta2 = NULL,
                               ncbiSummary = NULL,
                               prefix = NULL,
                               externalSoftwarePaths = NULL){
  
  parseNcbiSummaryR = NULL
  if(is.null(refSeqFasta2)){
    if(is.null(refSeqFasta1)){
      parseNcbiSummaryR = paste(c("Rscript", parseNcbiGeneR, 
                                  "-ncbi_summary", ncbiSummary, "-prefix", prefix, "-externalSoftware", externalSoftwarePaths), 
                                collapse=" ")
      
    }else{
      message(paste0("refseq_fasta1: ", refSeqFasta1))
      parseNcbiSummaryR = paste(c("Rscript", parseNcbiGeneR, "-refseq_fasta1", refSeqFasta1, 
                                  "-ncbi_summary", ncbiSummary, "-prefix", prefix, "-externalSoftware", externalSoftwarePaths), 
                                collapse=" ")
    }
  }else{
    parseNcbiSummaryR = paste(c("Rscript", parseNcbiGeneR, "-refseq_fasta1", refSeqFasta1, "-refseq_fasta2", 
                                refSeqFasta2,  "-ncbi_summary", ncbiSummary, "-prefix", prefix , "-externalSoftware", externalSoftwarePaths),
                              collapse=" ")    
  }
  system(parseNcbiSummaryR)
  
  # This is not ideal as this requires knowing the suffix
  # and how the names are generated before hand. It would
  # be better if this is some returned value.
  
  databaseName1 = paste(c(formatDir(getwd()), prefix, "_blastdb"), collapse="")
  databaseName2 = NULL
  if(length(refSeqFasta2) == 0){
    ##TODO: check with Sarah about the naming for the second database
    databaseName2 = paste(c(formatDir(getwd()), prefix, "_2_blastdb"), collapse="") 
    
  }
  ncbiAnnotFile = paste(c(formatDir(getwd()), prefix, ".ncbigene_annotation.txt"), collapse="")
  
  return(list(databaseName1 = databaseName1, databaseName2 = databaseName2, ncbiAnnotFile = ncbiAnnotFile))
}


##############################################################################################################
## Run kmer annotation
## @kmerAnnotationR: The path of R script kmerAnnotation.R
## @kmerAnalysisOutput: A list object containing the paths of the output files created by kmerAnalysis.R
## @dataBasePaths: A list object containing the paths of databases created by parse_ncbi_gene.R
## @ntBlastDB: A pre-created database
## @signif: The number of top significant kmers
## @nproc: The number of processors
## @prefix: The prefix of the output files created
##############################################################################################################
runKmerAnnotation = function(kmerAnnotationR = NULL,
                             kmerAnalysisOutput = NULL,
                             db1 = NULL,
                             db2 = NULL,
                             ncbiAnnot = NULL,
                             signif = NULL,
                             nproc = NULL,
                             prefix = NULL,
                             externalSoftwarePaths = NULL){
  
  chisqStat = kmerAnalysisOutput$chisqStat
  kmerResults = kmerAnalysisOutput$kmer
  nPresentCtrl = kmerAnalysisOutput$nPresentCtrl
  nPresentCase = kmerAnalysisOutput$nPresentCase
  blastdb1 = db1
  blastdb2 = db2
  
  kmerAnnotationCommand = paste(c("Rscript", kmerAnnotationR, 
                                  "-chisq_results", chisqStat, 
                                  "-kmer_results", kmerResults, 
                                  "-present_ctrl", nPresentCtrl,
                                  "-present_case", nPresentCase, 
                                  "-blastdb1", blastdb1, 
                                  "-blastdb2", blastdb2,
                                  "-ncbi_db", ncbiAnnot, 
                                  "-signif", signif, 
                                  "-nproc", nproc, 
                                  "-prefix", prefix, 
                                  "-externalSoftware", externalSoftwarePaths), 
                                collapse=" ")
  
  system(kmerAnnotationCommand)
  message("The kmers have been successfully annotated.")
  
}


###############################################################################################################
## Run kmer analysis with LMM
## @kmerLmmR: The path of the R script kmerLMM.R
## @kmerAnalysisOutput: A list object of all file paths of the output files created by kmerAnalysis.R
## @signifLMM: Number of top kmers to consider for the analysis with LMM
## @relateMatrix: The relation matrix created by GEMMA
## @prefix: Prefix of the output files
## @phenotypeFile: The path of the file containing the phenotype
###############################################################################################################
runKmerLmm = function(kmerLmmR = NULL,
                      kmerAnalysisOutput = NULL,
                      signifLMM = NULL,
                      relateMatrix = NULL,
                      prefix = NULL,
                      phenotypeFile = NULL,
                      externalSoftwarePaths = NULL){
  
  chisqStat = kmerAnalysisOutput$chisqStat
  patternKey = kmerAnalysisOutput$patternKey
  patternIndex = kmerAnalysisOutput$patternIndex
  
  lmmPrefix = paste(c(prefix,"lmm"), collapse="_")
  
  kmerLmmCommand = paste(c("Rscript", kmerLmmR,
                           "-chisqStat", chisqStat,
                           "-patternKey", patternKey,
                           "-patternIndex", patternIndex,
                           "-signif", signifLMM,
                           "-relateMatrix", relateMatrix, 
                           "-phenotype", phenotypeFile, 
                           "-prefix", lmmPrefix,
                           "-externalSoftware", externalSoftwarePaths), 
                         collapse=" ")
  message(kmerLmmCommand)
  system(kmerLmmCommand)  
  
  kmersUsed = paste(lmmPrefix, "_kmers_used.txt", sep="")
  lmmAllKmers = paste(lmmPrefix, "_LMM_allkmers_out.txt", sep="")
  kmerLmmOuput = list(kmersUsed = kmersUsed, lmmAllKmers = lmmAllKmers)
  
  return(kmerLmmOuput)
  
}


##############################################################################################################
## Run kmer annotation after the analysis with LMM
## @lmmKmerAnnotationR: The path of the R script LMM_kmer_annotation.R
## @kmerAnalysisOutput: A list object of all file paths of the output files created by kmerAnalysis.R 
## @dataBasePaths: A list object containing the paths of databases created by parse_ncbi_gene.R
## @ntBlastDB: A pre-created database
## @signif: The number of top significant kmers
## @nproc: The number of processors
## @prefix: The prefix of the output files created
## @lmmKmers: The path of the file containing the index of kmers used for LMM analysis (ends kmers_used.txt)
## @lmmOutput: The path of the file containing the Kmer LMM output for all kmers (ends allkmers_out.txt) 
##############################################################################################################
runLmmKmerAnnotation = function(lmmKmerAnnotationR = NULL,
                                kmerAnalysisOutput = NULL,
                                db1 = NULL,
                                db2 = NULL,
                                ncbiAnnot = NULL,
                                signif = NULL,
                                nproc = NULL,
                                prefix = NULL,
                                lmmKmers = NULL,
                                lmmOutput = NULL,
                                externalSoftwarePaths = NULL){
  
  chisqStat = kmerAnalysisOutput$chisqStat
  kmerResults = kmerAnalysisOutput$kmer
  nPresentCtrl = kmerAnalysisOutput$nPresentCtrl
  nPresentCase = kmerAnalysisOutput$nPresentCase
  
  lmmPrefix = paste(c(prefix,"lmm"), collapse="_")
  
  lmmKmerAnnotationCommand = paste(c("Rscript", lmmKmerAnnotationR,
                                     "-chisq_results", chisqStat, "-kmer_results", kmerResults, 
                                     "-present_ctrl", nPresentCtrl, "-present_case", nPresentCase, 
                                     "-blastdb1", db1, "-blastdb2", db2, "-ncbi_db", ncbiAnnot, 
                                     "-signif", signif, "-nproc", nproc, "-prefix", lmmPrefix, 
                                     "-LMM_kmers", lmmKmers, "-LMM_output", lmmOutput, "-externalSoftware", externalSoftwarePaths),
                                   collapse = " ")
  system(lmmKmerAnnotationCommand)
  
  
}


##### Read in the arguments from the command line #####
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


##### Get inputs from the command line #####
message("Arg length:")
message(length(args))
inputs = getCommandLineInputMatrix(args = args)
message(nrow(inputs))
message(ncol(inputs))

##### Input arguments required in general #####
dataFilePath = extractInputArgument(argName = "dataFile", commandLineInputs = inputs)
srcDir = formatDir(extractInputArgument(argName = "srcDir", commandLineInputs = inputs, default = getwd()))
prefix = extractInputArgument(argName = "prefix", commandLineInputs = inputs)
runLMM = as.logical(extractInputArgument(argName = "runLMM", commandLineInputs = inputs, default = TRUE))
geneticFileFormat = extractInputArgument(argName = "genFileFormat", commandLineInputs = inputs)
minCov = extractInputArgument(argName = "minCov", commandLineInputs = inputs, default = 5)
#toCreateKmerFiles = extractInputArgument(argName = "createKmerFiles", commandLineInputs = inputs, default = TRUE)
externalSoftwarePaths = extractInputArgument(argName="externalSoftware", commandLineInputs = inputs)
createDB = extractInputArgument(argName = "createDB", commandLineInputs = inputs, default = TRUE)
signif = extractInputArgument(argName = "signif", commandLineInputs = inputs, default = 10000)


##### Input arguments required for creating bam files #####
kmerFileDir = formatDir(extractInputArgument(argName = "kmerFileDir", commandLineInputs = inputs, default = getwd()))
qualityScore = extractInputArgument(argName = "qualityScore", commandLineInputs = inputs, default = "phred33")
nproc = as.numeric(extractInputArgument(argName = "nproc", commandLineInputs = inputs, default = 1))


###### Input specific to parse_ncbi_gene.R ######
refSeqFasta1 = extractInputArgument(argName = "refSeqFasta1", commandLineInputs = inputs, canBeNULL = TRUE)
refSeqFasta2 = extractInputArgument(argName = "refSeqFasta2", commandLineInputs = inputs, canBeNULL = TRUE)
ncbiSummary = extractInputArgument(argName = "ncbiSummary", commandLineInputs = inputs, canBeNULL = TRUE)


###### Input specific to kmerAnnotation.R ######
db1 = extractInputArgument(argName = "db1", commandLineInputs = inputs, canBeNULL = TRUE)
db2 = extractInputArgument(argName = "db2", commandLineInputs = inputs, canBeNULL = TRUE)
ncbiAnnot = extractInputArgument(argName = "ncbiAnnot", commandLineInputs = inputs, canBeNULL = TRUE)

###### Input specific to kmerLMM.R ######
relateMatrix = extractInputArgument(argName = "relateMatrix", commandLineInputs = inputs)
signifLMM = as.numeric(extractInputArgument(argName = "signifLMM", commandLineInputs = inputs, default = 100000))

message("Input arguments has been read.")

## Get the path of external software used
externalSoftwarePaths.df = read.table(file=externalSoftwarePaths, header=T, as.is = T)

parallelPath = getSoftwarePath("parallel", externalSoftwarePaths.df)


###### R scripts required ######
cat(paste0("Script directory: ",srcDir))
## Get the paths of the R scripts
doAndSortDskR = paste(srcDir, "doAndSortDsk.R", sep="")
createKmersFromFastaR = paste(srcDir, "createKmersFromFasta.R", sep="")
kmerAnalysisR = paste(srcDir, "kmerAnalysis.R", sep="")
parseNcbiGeneR = paste(srcDir, "parse_ncbi_gene.R", sep="")
kmerAnnotationR = paste(srcDir, "kmerAnnotation.R", sep="")
kmerLmmR = paste(srcDir, "kmerLMM.R", sep="")
lmmKmerAnnotationR = paste(srcDir, "LMM_KmerAnnotation.R", sep="")
## Check the paths of the R scripts
checkExistence(c(doAndSortDskR, kmerAnalysisR, parseNcbiGeneR, kmerAnnotationR, kmerLmmR, lmmKmerAnnotationR))


##### Preliminary checks ##### 
checkInputRelatedToDB(createDB = createDB, 
                      refSeqFasta1 = refSeqFasta1, refSeqFasta2 = refSeqFasta2, ncbiSummary = ncbiSummary, 
                      db1 = db1, db2 = db2, ncbiAnnot = ncbiAnnot)


data.df = read.table(file=dataFilePath, header=T, as.is=T)
phenotype = data.df$phenotype
toCreateKmerFiles = TRUE
geneticFileFormat = toupper(geneticFileFormat)
if(geneticFileFormat == "KMER"){
	toCreateKmerFiles = FALSE
}else if(!any(geneticFileFormat == c("BAM", "FASTA"))){
	stop("The genetic file format must be BAM, FASTA or Kmer.")
}


## Check if kmer files need to created
kmerPhenoPath = NULL
if(toCreateKmerFiles){
  
  bam_files = data.df$filePath
  comid = data.df$id
  
  ## Check that all the bam files exists
  if(geneticFileFormat == "BAM"){
  	checkExistence(bam_files)
  	kmerPhenoPath = createKmerFilesFromBAM(bam_files = bam_files, comid = comid, phenotype = phenotype,
    	 	                             nproc = nproc, kmerFileDir = kmerFileDir, prefix = prefix, 
        		                          externalSoftwarePaths = externalSoftwarePaths)
  }else{
  	kmerPhenoPath = createKmerFilesFromFASTA(createKmersFromFastaR = createKmersFromFastaR, dataFilePath = dataFilePath, 
  							kmerFileDir =  kmerFileDir, extSoft = externalSoftwarePaths)
  	
  }
  
}else{
  kmerPhenoPath = dataFilePath
  
}


##### Run kmer analysis #####
kmerAnalysisOutput = runKmerAnalysis(kmerAnalysisR = kmerAnalysisR, kmerPhenoPath = kmerPhenoPath, minCov = minCov,
                                         prefix = prefix, removeKmerTxt = FALSE, externalSoftwarePaths = externalSoftwarePaths)

##### Create BLAST databases #####
if(createDB){
  
  dataBasePaths = runParseNcbiSummary(parseNcbiGeneR = parseNcbiGeneR, refSeqFasta1 = refSeqFasta1, refSeqFasta2 = refSeqFasta2, 
                                      ncbiSummary = ncbiSummary, prefix = prefix, externalSoftwarePaths = externalSoftwarePaths)
  if(is.null(db1)){
    db1 = dataBasePaths$databaseName1
  }
  
  if(is.null(db2)){
    db2 = dataBasePaths$databaseName2
  }
  ncbiAnnot = dataBasePaths$ncbiAnnotFile
  
}

message(paste0("db1: ", db1))
message(paste0("db2: ", db2))
message(paste0("ncbiAnnot: ", ncbiAnnot))

runKmerAnnotation(kmerAnnotationR = kmerAnnotationR, kmerAnalysisOutput = kmerAnalysisOutput, 
                  db1 = db1, db2 = db2,  ncbiAnnot = ncbiAnnot,
                  signif = signif, nproc = nproc, prefix = prefix, externalSoftwarePaths = externalSoftwarePaths)



##### Run kmer analysis with LMM #####
if(runLMM){
  message(paste("Perform the kmer analysis with LMM:", runLMM, sep=" "))
  phenotypeFile = paste(c(formatDir(getwd()), prefix,"_phenotype.txt"), collapse="")
  message(paste(names(data.df)))
  write(data.df$phenotype, file=phenotypeFile, ncolumns = 1, append = FALSE, sep = "\n")

  kmerLmmOutput = runKmerLmm(kmerLmmR = kmerLmmR,
                             kmerAnalysisOutput = kmerAnalysisOutput, signifLMM = signifLMM, relateMatrix = relateMatrix, 
                             prefix = prefix, phenotypeFile = phenotypeFile, externalSoftwarePaths = externalSoftwarePaths)
  
 
  ##### Annotate kmers for analysis with LMM #####
  message("Perform kmer annotation for the analysis with LMM.")
  runLmmKmerAnnotation(lmmKmerAnnotationR = lmmKmerAnnotationR,
                       kmerAnalysisOutput = kmerAnalysisOutput,
                       db1 = db1, db2 = db2, ncbiAnnot = ncbiAnnot,
                       signif = signif, nproc = nproc, prefix = prefix, 
                       lmmKmers = kmerLmmOutput$kmersUsed, lmmOutput = kmerLmmOutput$lmmAllKmers,
                       externalSoftwarePaths = externalSoftwarePaths)
  
}

message("Kmer-based GWAS analysis completed successfully.")


