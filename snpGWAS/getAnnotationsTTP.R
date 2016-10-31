# Pull out extended annotation for a list of biallelic SNPs and genbank file
library("genoPlotR")
help = paste(
	"get_annotations.Rscript.R Pull out extended annotation for a list of biallelic SNPs and genbank file",
	"Daniel Wilson (2014)",
	"Chieh-Hsi Wu (2014)",
	"",
"Usage: Rscript get_annotations.Rscript.R snpinfo.txt reference.fasta ref_annnotation.gbk",
sep="\n")



#Genetic code:
geneticCode = list(
"TTT"="Phe","TTC"="Phe","TTA"="Leu","TTG"="Leu",
"TCT"="Ser","TCC"="Ser","TCA"="Ser","TCG"="Ser",
"TAT"="Tyr","TAC"="Tyr","TAA"="STO","TAG"="STO",
"TGT"="Cys","TGC"="Cys","TGA"="STO","TGG"="Trp",
"CTT"="Leu","CTC"="Leu","CTA"="Leu","CTG"="Leu",
"CCT"="Pro","CCC"="Pro","CCA"="Pro","CCG"="Pro",
"CAT"="His","CAC"="His","CAA"="Gln","CAG"="Gln",
"CGT"="Arg","CGC"="Arg","CGA"="Arg","CGG"="Arg",
"ATT"="Ile","ATC"="Ile","ATA"="Ile","ATG"="Met",
"ACT"="Thr","ACC"="Thr","ACA"="Thr","ACG"="Thr",
"AAT"="Asn","AAC"="Asn","AAA"="Lys","AAG"="Lys",
"AGT"="Ser","AGC"="Ser","AGA"="Arg","AGG"="Arg",
"GTT"="Val","GTC"="Val","GTA"="Val","GTG"="Val",
"GCT"="Ala","GCC"="Ala","GCA"="Ala","GCG"="Ala",
"GAT"="Asp","GAC"="Asp","GAA"="Glu","GAG"="Glu",
"GGT"="Gly","GGC"="Gly","GGA"="Gly","GGG"="Gly")



mutation_type_summary<-function(
	type1="-",
	type2="-",
	type3="-",
	type4="-",
	refCodon="-",
	nonRefCodon1="-",
	nonRefCodon2="-",
	nonRefCodon3="-",
	nonRefCodon4="-",
	frame="-",
	codonPosition="-",
	refAA="-",
	nonRefAA1="-",
	nonRefAA2="-",
	nonRefAA3="-",
	nonRefAA4="-"){
		
	ref=c(
		"Type1" = type1,
		"Type2" = type2,
		"Type3" = type3,
		"Type4" = type4,  
		"Refcodon" = refCodon, 
		"Nonrefcodon1" = nonRefCodon1, 
		"Nonrefcodon2" = nonRefCodon2, 
		"Nonrefcodon3" = nonRefCodon3,
		"Nonrefcodon4" = nonRefCodon4,  
		"Frame" = frame, 
		"Codonposition" = codonPosition, 
		"Refaa" = refAA, 
		"Nonrefaa1" = nonRefAA1,
		"Nonrefaa2" = nonRefAA2,
		"Nonrefaa3" = nonRefAA3,
		"Nonrefaa4" = nonRefAA4
	)
	
	return(ref)

}

getGeneticCode<-function(codonVec = NULL){
	
	aa<-as.vector(unlist(geneticCode[paste(codonVec, collapse="")]))
	if(is.null(aa)){
		aa<-NA
	}
	return(aa)
	
}

getMutationType<-function(wtaa = NULL, mtaa = NULL){
	type = NA
			
	if(is.na(wtaa) | is.na(mtaa)){			
		type = NA						
	}else if(wtaa==mtaa) {
		type = "Synonymous"
	} else if(wtaa=="STO") {
		type = "Read-through"
	} else if(mtaa=="STO") {
		type = "Nonsense"
	} else {
		type ="Non-synonymous"			
	}
	
	return(type)	
}

codonVecToString<-function(codonVec = NULL){
	codonStr<-paste(codonVec, collapse="")	
	return(codonStr)
}

# Crude estimate of a SNPs effect on protein function. 
# NB: snpinfo is a local variable so switching reference 
# and non-reference alleles has no effect. 
# Change to global variable if this behaviour is needed.

revcompl = c("A"="T", "C"="G", "G"="C", "T"="A")

mutation_type = function(i, ref_fa, snpinfo, ref_gbk, assert_Allele0=FALSE, logFilePath) {
		
	pos = as.numeric(snpinfo$Position[i])
	if(is.na(pos)) {
		return(mutation_type_summary(type1="Unknown", type2="Unknown", type3="Unknown", type4 = "Unknown"))
	}
	
	wh = which((ref_gbk$feature=="CDS" & ref_gbk$start<=pos & ref_gbk$end>=pos))
	if(length(wh)==0) {
		return(mutation_type_summary(
			type1="Intergenic", type2="Intergenic", type3="Intergenic", type4="Intergenic"))
	}
	
	wh = wh[which.max(ref_gbk$length[wh])]
	wt = NULL
	mt = NULL
	codonPosition = NULL
	frame = NULL
	if(ref_gbk$strand[wh]==1) {
		
		##Determine reading frame
		frame = (pos-ref_gbk$start[wh]) %% 3
		##Get wild type codon from the reference
		wt = ref_fa[pos-frame+0:2]
		
		## Check - to ensure major allele is also the consensus allele.
		## This is called if the reference base does not match with the major allele.
		if(wt[frame+1] != snpinfo$Allele0[i]) {
			
			if(assert_Allele0) {
				
				cat(paste("Warning 1: unexpected reference base at ",pos), 
					file = logFilePath, append=TRUE, sep = "\n")
					
			} else {
				
				cat(paste("Warning 2: unexpected reference base at ",pos), 
					file = logFilePath, append=TRUE, sep = "\n")
				
				wtFreq<-snpinfo[i, wt[frame+1]]	
				if(is.na(wtFreq)){
					
					cat(paste("Warning 3: ref at",pos, "is none of A, C, G or T.", sep=" ", collapse=""), 
						file = logFilePath, append = TRUE, sep = "\n")
						
				}else if(wtFreq == 0){
					
					cat(paste("Warning 4: ref base at", pos, "does not appear in the data."), 
						file = logFilePath, append = TRUE, sep = "\n")					
					
				}
				
			}
		}
		
		## The mutatant type is the minor allele
		mt<-matrix(rep(wt,4), ncol=3, byrow=T)
		allele.col.names<-c('Allele0','Allele1','Allele2','Allele3')
		
		for(j in 1:length(allele.col.names)){	
			allele = snpinfo[i, allele.col.names[j]]
			alleleFreq = snpinfo[i, allele]
			
			if(alleleFreq > 0){
				mt[j,frame+1] = allele
			}else{
				mt[j,]<-c('-','-','-')
			}
		}	
		
		codonPosition = ((pos-ref_gbk$start[wh]) %/% 3) + 1
		
		
	} else {
		cat(paste("The reverse compliment is interpreted at base idnex", pos, "."), 
						file = logFilePath, append = TRUE, sep = "\n")	
		
		frame = (ref_gbk$end[wh]-pos) %% 3
		wt = revcompl[ref_fa[pos+frame-(0:2)]]
		
		# Check
		if(wt[frame+1]!=revcompl[snpinfo$Allele0[i]]) {
			if(assert_Allele0) {
				
				cat(paste("Warning 5: unexpected reference base at ",pos), 
					file=logFilePath, append=TRUE, sep = "\n")
					
			} else {
				
				cat(paste("Warning 6: unexpected reference base at ",pos), 
					file = logFilePath, append=TRUE, sep = "\n")
				#convert the wild type back to the raw/reversed sequence
				wtFreq<-snpinfo[i, revcompl[wt[frame+1]]]	
				
				if(is.na(wtFreq)){
					
					cat(paste("Warning 7: ref at", pos, "is none of A, C, G or T.", sep=" ", collapse=""), 
						file = logFilePath, append = TRUE, sep = "\n")
						
				}else if(wtFreq == 0){
					
					cat(paste("Warning 8: ref base at", pos, "does not appear in the data."), 
						file = logFilePath, append = TRUE, sep = "\n")					
					
				}			
				
			}
		}
		
		
		## The mutatant type is the minor allele
		mt<-matrix(rep(wt,4), ncol=3, byrow=T)
		allele.col.names<-c('Allele0','Allele1','Allele2','Allele3')
		
		for(j in 1:length(allele.col.names)){	
			allele = snpinfo[i, allele.col.names[j]]
			alleleFreq = snpinfo[i, allele]
			
			if(alleleFreq > 0){
				mt[j,frame+1] = revcompl[allele]
			}else{
				mt[j,]<-c('-','-','-')
			}
		}
		
		codonPosition = ((ref_gbk$end[wh]-pos) %/% 3)+1
		
		
	}
	
	if(is.null(wt)){
		stop("Error: Wild type is not instantiated!")		
	}
	
	if(is.null(mt)){
		stop("Error: Mutant type is not instantiated!")
	}
	
	if(is.null(frame)){
		stop("Error: Frame is not instantiated!")
		
	}
	
	if(is.null(codonPosition)){
		stop("Error: Codon position is not instantiated!")
		
	}
	
	refCodon = codonVecToString(wt)
	nonRefCodonVec = apply(mt, 1, codonVecToString)	
		
	
		
	wtaa<-getGeneticCode(wt)
	mtaaVec = apply(mt, 1, getGeneticCode)		
	mutTypes = sapply(mtaaVec, getMutationType, wtaa=wtaa)
	names(mutTypes)<-NULL
	names(nonRefCodonVec)<-NULL
	names(mtaaVec)<-NULL
		

	return(mutation_type_summary(
		type1=mutTypes[1],
		type2=mutTypes[2],
		type3=mutTypes[3],
		type4=mutTypes[4],
		refCodon=refCodon,
		nonRefCodon1=nonRefCodonVec[1],
		nonRefCodon2=nonRefCodonVec[2],
		nonRefCodon3=nonRefCodonVec[3],
		nonRefCodon4=nonRefCodonVec[4],
		frame=frame+1,
		codonPosition=codonPosition,
		refAA=wtaa,
		nonRefAA1=mtaaVec[1],
		nonRefAA2=mtaaVec[2],
		nonRefAA3=mtaaVec[3],
		nonRefAA4=mtaaVec[4])	
	)	
}


# Identify the genes each variant lands in or between
gbk.locate = function(gbk,position) {
	position = as.numeric(position)
	beg = which(gbk$start<=position); beg = beg[length(beg)]
	end = which(position<=gbk$end)[1]
	nomatch = length(beg)==0 | length(end)==0
	if(!nomatch) nomatch = nomatch | is.na(beg) | is.na(end)
	if(nomatch) {
		tp = gbk[1,]
		tp[1,1:ncol(tp)] = "-"
		return(tp)
	}
	as.data.frame(t(apply(gbk[beg:end,,drop=FALSE],2,paste0,collapse=":")))
}


# Assumes a fasta file representing a single genome, possibly split across contigs
read.fasta.ref = function(ref_file) {
	r = scan(ref_file,what=character(0),sep="\n")
	beg = substr(r,1,1)
	gd = beg!=">"
	rcat = paste(r[gd],collapse="")
	return(toupper(unlist(strsplit(rcat,""))))
}

getAnnotations<-function(
	snpinfo_file = NULL,
	reference_file = NULL,
	genbank_file = NULL, 
	prefix = NULL){
		
	snpinfo = read.delim(snpinfo_file, header=T, as.is=TRUE, sep="\t", quote="", comment="")
	ref = read.fasta.ref(reference_file)
	ref_gbk_pipe<-pipe(paste0("sed 's/\\/note=\"\\*matching_locus_tag: /\\/matching_locus_tag=\"/g' ", genbank_file))
	ref_gbk = read_dna_seg_from_file(ref_gbk_pipe, fileType="genbank", extra_fields="matching_locus_tag")

	#splitn=floor(nrow(snpinfo)/10)
	splitseq=seq(1, by=1000, to=nrow(snpinfo))
	splitseq[length(splitseq)+1]=(nrow(snpinfo)+1)

	header=c(## From snpinfo.txt
		names(snpinfo), 
		## From mutation_type
		"Type1", "Type2", "Type3", "Type4", 
		"Refcodon", "Nonrefcodon1", "Nonrefcodon2", "Nonrefcodon3", "Nonrefcodon4",
		"Frame", "Codonposition", 
		"Refaa", "Nonrefaa1", "Nonrefaa2", "Nonrefaa3", "Nonrefaa4",
		## From gbk.locate 
		"name", "start", "end", "strand", "length", "pid", 
		"gene", "synonym", "product", "proteinid", "feature", "gene_type", 
		"matching_locus_tag", "col", "lty", "lwd", "pch", "cex")
	
	logFilePath = paste0(prefix,"_ttp.annotation_logfile.txt")
	annotFileName = paste0(prefix,".ttp.annotation.txt")
	cat(paste(c(header),collapse="\t"), file=annotFileName, append=FALSE, sep="\n")
		
	for(j in 1:(length(splitseq)-1)){
		
		snpinfo2 = snpinfo[splitseq[j]:(splitseq[(j+1)]-1),]
		annot = list()
		
		for(i in 1:nrow(snpinfo2)) {
			
			annot1 = t(mutation_type(i, ref, snpinfo2, ref_gbk, assert_Allele0=FALSE, logFilePath))
			annot2 = gbk.locate(ref_gbk, snpinfo2$Position[i])
			new.annot<-cbind(snpinfo2[i,], annot1, annot2[1,])
			annot.names<-names(annot)
			new.annot.names<-names(new.annot)
			
			
			
			annot = rbind(annot,new.annot)
		}
		cat("Done",j*1000,sep="\n")
		
		
		write.table(annot,file=annotFileName,row=F,col=F,quote=F,sep="\t",append = TRUE)
		cat(paste0(1000*j, " snps annotated: ",Sys.time()), 
			file = logFilePath, append=(j > 1), sep="\n")
	}
	
	close(ref_gbk_pipe)
	
	return(annotFileName)
	
	
}


# Read options from command line
args = commandArgs(trailingOnly = TRUE)
if(!(length(args) ==4 | length(args) == 0)) {
	cat(help,sep="\n")
	stop("\nIncorrect usage\n")
}

if(length(args) == 4){
	snpinfo_file = args[1]
	#snpinfo_file = "/Users/jessiewu/Documents/gwas/ecoli/snpAnalysis/multiAllelic/ecoli245.snpinfo.1000.txt"
	
	reference_file = args[2]
	#reference_file = "/Users/jessiewu/Documents/gwas/ecoli/snpAnalysis/R00000042.fa"
	
	genbank_file = args[3]
	#genbank_file = "/Users/jessiewu/Documents/gwas/ecoli/snpAnalysis/R00000042.gbk"
	
	prefix=args[4]
	#prefix = "ecoli245"


	getAnnotations(
		snpinfo_file = snpinfo_file, 
		reference_file = reference_file, 
		genbank_file = genbank_file, 
		prefix = prefix
	)
}
