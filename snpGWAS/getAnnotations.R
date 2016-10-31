# Pull out extended annotation for a list of biallelic SNPs and genbank file
library("genoPlotR")
help = paste(
"get_annotations.Rscript.R Pull out extended annotation for a list of biallelic SNPs and genbank file",
"Daniel Wilson (2014)",
"",
"Usage: Rscript get_annotations.Rscript.R bipinfo.txt reference.fasta ref_annnotation.gbk",
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


# Crude estimate of a SNPs effect on protein function
revcompl = c("A"="T","C"="G","G"="C","T"="A")

# Crude estimate of a SNPs effect on protein function. NB: bipinfo is a local variable so switching reference and non-reference alleles has no effect. Change to global variable if this behaviour is needed.
mutation_type = function(i,ref_fa,bipinfo,ref_gbk,assert_Allele0=FALSE,logFilePath) {
	ret = c("Type"="-","Refcodon"="-","Nonrefcodon"="-","Frame"="-","Codonposition"="-","Refaa"="","Nonrefaa"="")
	revcompl = c("A"="T","C"="G","G"="C","T"="A")
	pos = as.numeric(bipinfo$Position[i])
	if(is.na(pos)) {
		ret[1] = "Unknown"
		return(ret)
	}
	wh = which((ref_gbk$feature=="CDS" & ref_gbk$start<=pos & ref_gbk$end>=pos))
	if(length(wh)==0) {
		ret[1] = "Intergenic"
		return(ret)
	}
	wh = wh[which.max(ref_gbk$length[wh])]
	if(ref_gbk$strand[wh]==1) {
		frame = (pos-ref_gbk$start[wh]) %% 3
		wt = ref_fa[pos-frame+0:2]
		# Check
		if(wt[frame+1]!=bipinfo$Allele0[i]) {
			if(assert_Allele0) {
				
				cat(paste("Warning 1: unexpected reference base at ",pos), 
					file=logFilePath, append=TRUE, sep = "\n")
			} else {
				cat(paste("Warning 2: unexpected reference base at ",pos), 
					file=logFilePath, append=TRUE, sep = "\n")
				if(wt[frame+1]!=bipinfo$Allele1[i]) {
					cat(paste("Warning 3: neither reference nor non-reference base match at ",pos), 
						file=logFilePath, append=TRUE, sep = "\n")
				} else {
					tmp = bipinfo$Allele1[i]
					bipinfo$Allele1[i] = bipinfo$Allele0[i]
					bipinfo$Allele0[i] = tmp
				}
			}
		}
		mt = wt; mt[frame+1] = bipinfo$Allele1[i]
		wtaa = as.vector(unlist(geneticCode[paste(wt,collapse="")]))
		mtaa = as.vector(unlist(geneticCode[paste(mt,collapse="")]))
		if(length(wtaa) == 0 | length(mtaa) == 0){
			ret = c("Type" = NA,
				"Refcodon"=paste(wt,collapse=""),
				"Nonrefcodon"=paste(mt,collapse=""),
				"Frame"=frame+1,
				"Codonposition"=((pos-ref_gbk$start[wh]) %/% 3)+1,
				"Refaa"=NA,
				"Nonrefaa"=NA
				)
			return (ret)
			
		}
		
		if(length(wtaa) == 0 ){
			message(i)
			message(paste(wt,collapse=""))			
		}
		if(length(mtaa) == 0){
			message(i)
			message(paste(mt,collapse=""))			
		}
		if(wtaa==mtaa) {
			
			ret = c("Type"="Synonymous",
				"Refcodon"=paste(wt,collapse=""),
				"Nonrefcodon"=paste(mt,collapse=""),
				"Frame"=frame+1,
				"Codonposition"=((pos-ref_gbk$start[wh]) %/% 3)+1,
				"Refaa"=wtaa,
				"Nonrefaa"=mtaa
				)
			return(ret)
			
		} else if(wtaa=="STO") {
			
			ret = c("Type"="Read-through",
				"Refcodon"=paste(wt,collapse=""),
				"Nonrefcodon"=paste(mt,collapse=""),
				"Frame"=frame+1,
				"Codonposition"=((pos-ref_gbk$start[wh]) %/% 3)+1,
				"Refaa"=wtaa,
				"Nonrefaa"=mtaa
				)
			return(ret)
			
		} else if(mtaa=="STO") {
			
			ret = c("Type"="Nonsense",
				"Refcodon"=paste(wt,collapse=""),
				"Nonrefcodon"=paste(mt,collapse=""),
				"Frame"=frame+1,
				"Codonposition"=((pos-ref_gbk$start[wh]) %/% 3)+1,
				"Refaa"=wtaa,
				"Nonrefaa"=mtaa
				)
			return(ret)
			
		} else {
			
			ret = c("Type"="Non-synonymous",
			"Refcodon"=paste(wt,collapse=""),
			"Nonrefcodon"=paste(mt,collapse=""),
			"Frame"=frame+1,
			"Codonposition"=((pos-ref_gbk$start[wh]) %/% 3)+1,
			"Refaa"=wtaa,
			"Nonrefaa"=mtaa)
			return(ret)
			
		}
	} else {
		frame = (ref_gbk$end[wh]-pos) %% 3
		wt = revcompl[ref_fa[pos+frame-(0:2)]]
		# Check
		if(wt[frame+1]!=revcompl[bipinfo$Allele0[i]]) {
			if(assert_Allele0) {
				cat(paste("Warning 4: unexpected reference base at ",pos), 
					file=logFilePath, append=TRUE, sep = "\n")
			} else {
				cat(paste("Warning 5: unexpected reference base at ",pos), 
					file=logFilePath, append=TRUE, sep = "\n")
				if(wt[frame+1]!=revcompl[bipinfo$Allele1[i]]) {
					cat(paste("Neither reference nor non-reference base match at ", pos),
					 file=logFilePath, append=TRUE, sep = "\n")
				} else {
					tmp = bipinfo$Allele1[i]
					bipinfo$Allele1[i] = bipinfo$Allele0[i]
					bipinfo$Allele0[i] = tmp
				}
			}
		}
		mt = wt; mt[frame+1] = revcompl[bipinfo$Allele1[i]]
		wtaa = as.vector(unlist(geneticCode[paste(wt,collapse="")]))
		mtaa = as.vector(unlist(geneticCode[paste(mt,collapse="")]))
		if(length(wtaa) == 0 | length(mtaa) == 0){
			ret = c("Type" = NA,
				"Refcodon"=paste(wt,collapse=""),
				"Nonrefcodon"=paste(mt,collapse=""),
				"Frame"=frame+1,
				"Codonposition"=((pos-ref_gbk$start[wh]) %/% 3)+1,
				"Refaa"=NA,
				"Nonrefaa"=NA
				)
			return (ret)
			
		}
		if(wtaa==mtaa) {
			
			ret = c("Type"="Synonymous",
				"Refcodon"=paste(wt,collapse=""),
				"Nonrefcodon"=paste(mt,collapse=""),
				"Frame"=frame+1,
				"Codonposition"=((ref_gbk$end[wh]-pos) %/% 3)+1,
				"Refaa"=wtaa,
				"Nonrefaa"=mtaa
				)
			return(ret)
			
		} else if(wtaa=="STO") {
			ret = c("Type"="Read-through",
					"Refcodon"=paste(wt,collapse=""),
					"Nonrefcodon"=paste(mt,collapse=""),
					"Frame"=frame+1,
					"Codonposition"=((ref_gbk$end[wh]-pos) %/% 3)+1,
					"Refaa"=wtaa,
					"Nonrefaa"=mtaa
					)
			return(ret)
		} else if(mtaa=="STO") {
			ret = c("Type"="Nonsense",
					"Refcodon"=paste(wt,collapse=""),
					"Nonrefcodon"=paste(mt,collapse=""),
					"Frame"=frame+1,
					"Codonposition"=((ref_gbk$end[wh]-pos) %/% 3)+1,
					"Refaa"=wtaa,
					"Nonrefaa"=mtaa
					)
			return(ret)
		} else {
			ret = c("Type"="Non-synonymous",
					"Refcodon"=paste(wt,collapse=""),
					"Nonrefcodon"=paste(mt,collapse=""),
					"Frame"=frame+1,
					"Codonposition"=((ref_gbk$end[wh]-pos) %/% 3)+1,
					"Refaa"=wtaa,
					"Nonrefaa"=mtaa
					)
			return(ret)
		}
	}
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
	bipinfo_file = NULL,
	reference_file = NULL,
	genbank_file = NULL, 
	prefix = NULL){
		
	bipinfo = read.delim(bipinfo_file, header = TRUE, as.is=TRUE,sep="\t",quote="",comment="")
	ref = read.fasta.ref(reference_file)
	#genbank_file = "/Volumes/UDISK/svn/gwas_methods/release/gwasSourceCodes/example/R00000022.gbk"
	#genbank_file = "/Volumes/UDISK/svn/gwas_methods/release/gwasSourceCodes/example/sim/simRef.gbk"
	message(genbank_file)
	ref_gbk = read_dna_seg_from_file(pipe(paste0("sed 's/\\/note=\"\\*matching_locus_tag: /\\/matching_locus_tag=\"/g' ", genbank_file)), fileType="genbank", extra_fields="matching_locus_tag")

	#splitn=floor(nrow(bipinfo)/10)
	splitseq=seq(1,by=1000,to=nrow(bipinfo))
	splitseq[length(splitseq)+1]=(nrow(bipinfo)+1)

	header=c(names(bipinfo), "Type", "Refcodon", "Nonrefcodon", "Frame", "Codonposition", "Refaa", "Nonrefaa", "name", "start", "end", "strand", "length", "pid", "gene", "synonym", "product", "proteinid", "feature", "gene_type", "matching_locus_tag", "col", "lty", "lwd", "pch", "cex")
	cat(paste(c(header),collapse="\t"),file=paste0(prefix,".bip.annotation.txt"),append=FALSE,sep="\n")
	logFilePath = paste0(prefix,"_bipannotation_logfile.txt")
	annotFileName = paste0(prefix,".bip.annotation.txt")
	for(j in 1:(length(splitseq)-1)){
		
		bipinfo2 = bipinfo[splitseq[j]:(splitseq[(j+1)]-1),]
		annot = list()
		
		for(i in 1:nrow(bipinfo2)) {
			
			annot1 = t(mutation_type(i,ref,bipinfo2,ref_gbk,assert_Allele0=FALSE,logFilePath))
			annot2 = gbk.locate(ref_gbk,bipinfo2$Position[i])
			annot = rbind(annot,cbind(bipinfo2[i,],annot1,annot2[1,]))
			#if(nrow(annot2)>1){
			#	for(k in 2:nrow(annot2)){
			#		annot = rbind(annot,
			#					cbind("Position"="-", "Allele0"="-", "Allele1"="-", annot1, annot2[k,])
			#				)
			#	}
			#}
		}
		cat("Done",j*1000,sep="\n")
		
		
		write.table(annot,file=annotFileName,row=F,col=F,quote=F,sep="\t",append=TRUE)
		cat(paste0(1000*j, " bips annotated: ",Sys.time()), 
			file = logFilePath, append=(j > 1), sep="\n")
	}
	
	return(annotFileName)
	
	
	
}


# Read options from command line
args = commandArgs(trailingOnly = TRUE)
if(!(length(args) ==4 | length(args) == 0)) {
	cat(help,sep="\n")
	stop("\nIncorrect usage\n")
}
if(length(args) ==4){
	bipinfo_file = args[1]
	#bipinfo_file = "/Users/jessiewu/Documents/gwas/ecoli/snpAnalysis/ecoli245.bipinfo.txt"
	reference_file = args[2]
	#reference_file = "/Users/jessiewu/Documents/gwas/ecoli/snpAnalysis/R00000042.fa"
	genbank_file = args[3]
	#genbank_file = "/Users/jessiewu/Documents/gwas/ecoli/snpAnalysis/R00000042.gbk"
	prefix=args[4]
	#prefix = "ecoli245"


	getAnnotations(
		bipinfo_file = bipinfo_file, 
		reference_file = reference_file, 
		genbank_file = genbank_file, 
		prefix = prefix
	)
}


#bipinfo = read.delim(bipinfo_file,as.is=TRUE,sep="\t",quote="",comment="")
#bipinfo2 = bipinfo[43001:44000,]
#ref = read.fasta.ref(reference_file)
#ref_gbk = read_dna_seg_from_file(pipe(paste0("sed 's/\\/note=\"\\*matching_locus_tag: /\\/matching_locus_tag=\"/g' ", genbank_file)), fileType="genbank", extra_fields="matching_locus_tag")
#i=711
#annot1 = t(mutation_type(711,ref,bipinfo2,ref_gbk,assert_Allele0=FALSE,logFilePath))
#			annot2 = gbk.locate(ref_gbk,bipinfo2$Position[i])
#cbind(bipinfo2[i,],annot1,annot2[1,])