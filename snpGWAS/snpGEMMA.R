runGEMMA<-function(
	prefix = NULL,
	gemmaPath = NULL,
	gemmaGENFormatFilePath = NULL,
	gemmaSNPFormatFilePath = NULL,
	filepathsPhenoTablePath = NULL){
	
	
	phenotype_file = paste0(prefix,"phenotype.txt")
	filepathsPhenoTable = read.table(file=filepathsPhenoTablePath,header=T)
	write.table(filepathsPhenoTable$phenotype,
		file = phenotype_file,
		quote = FALSE,
		row.names = FALSE,
		col.names = FALSE
	)
	
	
	## Get relatedness matrix
	system(paste0(gemmaPath, " -g ", gemmaGENFormatFilePath,
					" -p ", phenotype_file,
					" -a " , gemmaSNPFormatFilePath ,
					" -gk 1 -o ", prefix,"_gemma_relmatrixout -maf 0"
					)
			)
	
	## LMM
	system(paste0(gemmaPath, 
					" -g ", gemmaGENFormatFilePath ,
					" -p ", phenotype_file,
					" -a ", gemmaSNPFormatFilePath, 
					" -k ./output/",prefix,"_gemma_relmatrixout.cXX.txt -lmm 4 -o ",
					prefix,"_gemma_lmmout -maf 0"
				)
			)


	## Read output
	gem=read.table(paste0("./output/",prefix,"_gemma_lmmout.assoc.txt"), header=TRUE, as.is=TRUE)
	#SNPid=gem$rs
	position=gem$ps
	#pval=as.numeric(gem$p_score)
	#pvallik=gem$p_lrt
	#gem=cbind(SNPid,position,pvallik)
	pvalWald = as.numeric(gem$p_wald)
	pvalLRT = as.numeric(gem$p_lrt)
	pvalScore = as.numeric(gem$p_score)
	negLog10PValWald=-log10(pvalWald)
	negLog10PValLRT=-log10(pvalLRT)
	negLog10PValScore=-log10(pvalScore)
	#message("pval:")
	#message(paste(pval[1],pval[2]))
	#message(paste(log[1],log[2]))
	#gem=cbind(gem,log)
	
	#data5=subset(gem,gem[,4] >=signif)

	

	# manhattan.plot(
		# position = position, 
		# pvalue = negLog10PValWald,
		# file = paste0(prefix,"_SNP_gemma_manhattan_wald.png"),
		# plot.title = paste0(prefix," SNP GWAS gemma manhatten plot (Wald test)")
	# )
	
	# manhattan.plot(
		# position = position, 
		# pvalue = negLog10PValLRT,
		# file = paste0(prefix,"_SNP_gemma_manhattan_lrt.png"),
		# plot.title = paste0(prefix," SNP GWAS gemma manhatten plot (LRT)")
	# )
	
	# manhattan.plot(
		# position = position, 
		# pvalue = negLog10PValScore,
		# file = paste0(prefix,"_SNP_gemma_manhattan_score.png"),
		# plot.title = paste0(prefix," SNP GWAS gemma manhatten plot (z-score)")
	# )

	cat(paste0("GEMMA completed: ",Sys.time()),file=paste0(prefix,"_input_logfile.txt"),append=TRUE,sep="\n")
	
	gemmaAnalysisInfo<-list(SNPid=gem$rs,
		position=position,
		pvalWald = pvalWald,
		pvalLRT = pvalLRT,
		pvalScore = pvalScore,
		negLog10PValWald = negLog10PValWald,
		negLog10PValLRT = negLog10PValLRT,
		negLog10PValScore = negLog10PValScore
	)
	
	return(gemmaAnalysisInfo)
	
}

processGEMMA<-function(
	SNPid = NULL,
	position = NULL,
	pvalWald = NULL,
	pvalLRT = NULL,
	pvalScore = NULL,
	negLog10PValWald = NULL,
	negLog10PValLRT = NULL,
	negLog10PValScore = NULL,
	annotationPath = NULL,
	prefix = NULL,
	signif.threshold = NULL){

	annot=scan(annotationPath, what=character(0),sep="\n")
	annot.header=annot[1]
	annot=annot[-1]
		
	manhattan.plot(
		position = position, 
		pvalue = negLog10PValWald,
		file = paste0(prefix,"_bip_gemma_manhattan_wald.png"),
		plot.title = paste0(prefix," (biallelic) SNP GWAS gemma manhatten plot (Wald)")
	)
	
	manhattan.plot(
		position = position, 
		pvalue = negLog10PValLRT,
		file = paste0(prefix,"_bip_gemma_manhattan_lrt.png"),
		plot.title = paste0(prefix," (biallelic) SNP GWAS gemma manhatten plot (LRT)")
	)
	
	manhattan.plot(
		position = position, 
		pvalue = negLog10PValScore,
		file = paste0(prefix,"_bip_gemma_manhattan_score.png"),
		plot.title = paste0(prefix," (biallelic) SNP GWAS gemma manhatten plot (z-score)")
	)
	
	
	gemmaSummary.df<-data.frame(
		SNPid = SNPid,
		position = position,
		pvalWald = pvalWald,
		pvalLRT = pvalLRT,
		pvalScore = pvalScore,
		negLog10PValWald = negLog10PValWald,
		negLog10PValLRT = negLog10PValLRT,
		negLog10PValScore = negLog10PValScore,
		annot
	)
	
	threshold.set=subset(gemmaSummary.df, gemmaSummary.df$negLog10PValLRT >= signif.threshold)
	colnames(threshold.set)[ncol(threshold.set)] = annot.header
	negLog10POrder = order(threshold.set$negLog10PValLRT, decreasing=T)
	write.table(threshold.set[negLog10POrder,], 
		file=paste0(prefix,"_gemma_signifannot.txt"), quote=F, row=F, sep="\t")
	write.table(gemmaSummary.df[,-ncol(gemmaSummary.df)], 
		file=paste0(prefix,"_gemma_analysis_full.txt"), quote=F, row=F, sep="\t")

}

snpAnalysisGEMMA<-function(
	prefix = NULL,
	gemmaPath = NULL,
	gemmaGENFormatFilePath = NULL,
	gemmaSNPFormatFilePath = NULL,
	filepathsPhenoTablePath = NULL,
	annotationPath = NULL,
	signif.threshold = NULL){
		
	
	
		
	gemmaAnalysisInfo<-runGEMMA(
		prefix = prefix,
		gemmaPath = gemmaPath,
		gemmaGENFormatFilePath = gemmaGENFormatFilePath,
		gemmaSNPFormatFilePath = gemmaSNPFormatFilePath,
		filepathsPhenoTablePath = filepathsPhenoTablePath
	)
	
	
	
	
	processGEMMA(
		SNPid = gemmaAnalysisInfo$SNPid,
		position = gemmaAnalysisInfo$position,
		pvalWald = gemmaAnalysisInfo$pvalWald,
		pvalLRT = gemmaAnalysisInfo$pvalLRT,
		pvalScore = gemmaAnalysisInfo$pvalScore,
		negLog10PValWald = gemmaAnalysisInfo$negLog10PValWald,
		negLog10PValLRT = gemmaAnalysisInfo$negLog10PValLRT,
		negLog10PValScore = gemmaAnalysisInfo$negLog10PValScore,
		annotationPath = annotationPath,
		prefix = prefix,
		signif.threshold = signif.threshold
	)
	
	
		
}




args = commandArgs(trailingOnly = TRUE)
message(length(args))
if(!(length(args) == 8 | length(args) ==0)) {
	stop("\nIncorrect usage\n")
}

if(length(args) == 8){
	source(args[8])
	snpAnalysisGEMMA(
		prefix = args[1],
		gemmaPath = args[2],
		gemmaGENFormatFilePath = args[3],
		gemmaSNPFormatFilePath = args[4],
		filepathsPhenoTablePath = args[5],
		annotationPath = args[6],
		signif.threshold = args[7]
	)

}
