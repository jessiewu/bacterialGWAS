#' This function retrieves the binary data from GEMMA input files.
#' 
#' This function performs the genome-wide principal components.
#' @param gemmaGenFile The path to the GEMMA gen input file. It is a required input.
#' @param gemmaSnpFile The path to the GEMMA snp input file. It is a required input.
#' @param id A vector of the sample ids.
#' @param prefix Prefix of the output files.
#' @keywords GEMMA
#' @keywords SNPs
#' @return A list of information on the biallelic snps
#' @export
#' @examples
#' pheno <- read.table(file = pheno, header=T, as.is = T, sep="\t")
#' data <- getDataFromGemma(gemmaGenFile = gemmaGenFile, gemmaSnpFile = gemmaSnpFile, 
#'								id = pheno.df$id, prefix = prefix)
getDataFromGemma = function(prefix = NULL, gemmaGenFile = NULL, gemmaSnpFile = NULL, id = NULL){
	
	gemmaGen.df = read.table(file=gemmaGenFile, header=F, as.is=T)
	gemmaSNP.df = read.table(file=gemmaGenFile, header=F, as.is=T)
	bipCount = nrow(gemmaGen.df)
	dataCount = ncol(gemmaGen.df) - 3
	
	if(!is.null(id)){
		if(length(id) != dataCount){
			stop("The number of id's (", length(id), ") does not match up with the number of samples(", dataCount,".")
		}
	}
	
	
	m = matrix(nrow=2,ncol=bipCount)
	gen = gemmaGen.df[,-c(1:3)]
	m[2,] = rowSums(gemmaGen.df[,-c(1:3)])
	m[1,] = dataCount - m[2,]

	allele.id <- matrix(c(0, 1)[apply(m, 2, order, decreasing=TRUE)], nrow=2)

	# Output filenames
	# biallelic polymorphisms encoded -1 (missing) 0 (allele 0) 1 (allele 1)
	bip_outfile <- paste0(prefix, ".gemma.bip.patterns.txt");		
	
	# positional and allelic information for biallelic polymorphisms
	bipinfo_outfile <- paste0(prefix, ".gemma.bipinfo.txt");		
	
	# Allocate memory for bip and snp, because need to transform so cannot output on the fly
	bip <- matrix(NA, bipCount, dataCount)

	for(i in 1:ncol(gen)) {
		# Read the mapcall file
		fa <- gen[,i]
		bip[, i] <- -1
		bip[fa==allele.id[1, ], i] = 0
		bip[fa==allele.id[2, ], i] = 1
		
	}

	# Convert BIP and SNP patterns to factors to identify equivalencies
	bip.pat <- factor(apply(bip, 1, paste, collapse=""))
	# Record only unique patterns, and record the pattern equivalence in the bipinfo file
	bip.pat1 <- match(levels(bip.pat), bip.pat)
	
	# Output compacted bip and snp objects
	if(is.null(id)){
		id = paste("id", c(1:ncol(bip)), sep="")
	}
	colnames(bip) = id
	write.table(bip[bip.pat1, ], bip_outfile, row=FALSE, col=TRUE, quote=FALSE, sep="\t")
	
	# Output info files
	pos = as.numeric(gemmaSNP.df[,2])
	bipinfo <- data.frame("Position"= pos,
						 "Allele0"=allele.id[1,],
						 "Allele1"=allele.id[2,],
						 "0"=m[1,],
						 "1"=m[2,],
						 "Pattern"=as.numeric(bip.pat));

	write.table(bipinfo, bipinfo_outfile, row=FALSE, col=TRUE, quote=FALSE, sep="\t")
	ps.bips <- bipinfo$Position
	bippat <- sapply(1:max(bipinfo$Pattern), function(x)sum(bipinfo$Pattern==x))
	ipat <- bipinfo$Pattern
	rm(bipinfo)		
	
	return(list("XX" = bip[bip.pat1, ],
				"XX.tritetra" = NULL,
				"bippat" = bippat,
				"snppat" = NULL,
				"pattern" = ipat,
				"pattern.snps" = NULL,
				"ps" = ps.bips,
				"ps.snps" = NULL,
				"n.triallelic" = 0,
				"n.tetraallelic" = 0)

}