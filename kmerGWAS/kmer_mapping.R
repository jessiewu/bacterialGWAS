# Annotate kmers by mapping to a reference using bowtie2
help=paste("","kmer_mapping.R","Usage: Rscript kmer_mapping.R bowtie_output reference_genbank reference_fasta output_prefix mapping_quality","","bowtie_output: Output from running bowtie on a list of kmers","reference_genbank: The genbank file of the reference used for bowtie mapping","reference_fasta: The fasta file of the reference used for bowtie mapping","output_prefix: The prefix to all output files","mapping_quality: Which bowtie quality score cut off to use (Optional, default=0, annotates all kmers)",sep="\n")

library(genoPlotR)

annotation = function(pos){
	sapply(pos, function(i) which(genbank$start <= i & genbank$end >= i))
}

read_reference = function(ref_file) {
	r = scan(ref_file,what=character(0),sep="\n")
	rcat = paste(r[2:length(r)],collapse="")
	return(unlist(strsplit(rcat,"")))
}


args = commandArgs(trailingOnly = TRUE)
if(length(args!=0)){
if(args[1]=="-help" | args[1]=="-h"){
	cat(help,sep="\n")
	q("no")
}
}

if(length(args)<4 | length(args)==0) {
		    cat(help,sep="\n")
		    stop("\nIncorrect usage\n")
}


kmerfile=args[1]
genbank=args[2]
reference_fa=args[3]
prefix=args[4]
if(length(args)==5){
	quality=as.numeric(args[5])
} else {
	quality=0
}

keep_unmapped=0

genbank=read_dna_seg_from_genbank(genbank)

nlines=as.numeric(scan(pipe(paste0("samtools view -S ",kmerfile," -F 0 -q ",quality," | wc -l")),what=character(0)))

len_results=floor(nlines/10000)

reference_fa=read_reference(reference_fa)

ref_length=length(reference_fa)

a=c()
for(i in 1:len_results){
	a=scan(pipe(paste0("samtools view -S ",kmerfile," -F ",keep_unmapped," -q ",quality," | head -n ",i,"0000"," | tail -n 10000")),what=character(0),sep="\n",quote="")
	bowtie_out=matrix(NA,ncol=23,nrow=nlines)
	for(j in 1:length(a)){
		t=unlist(strsplit(a[j],"\t"))
		bowtie_out[j,1:length(t)]=t
	}
	bowtie_out[,1]=as.numeric(bowtie_out[,1])+1
	pos=as.numeric(bowtie_out[,4])
	gen=annotation(pos)
	genes=genbank[sapply(1:length(gen), function(x) ifelse(length(gen [[x]]) == 0, NA, gen[[x[1]]])),]
	bowtie_out=cbind(bowtie_out,genes)
	len_nomatch=length(which(bowtie_out[,2]=="4"))
	new_pos=seq(from=(ref_length+25000),by=100,length=len_nomatch)
	bowtie_out[,4]=as.character(bowtie_out[,4])
	bowtie_out[,4]=as.numeric(bowtie_out[,4])
	bowtie_out[which(bowtie_out[,2]=="4"),4]=new_pos
	write.table(bowtie_out,file=paste0(prefix,"_bowtie_annot.txt"),row=F,col=F,quote=F,sep="\t",append=(i>1))
	cat(paste0("Done ",i,"0000"),sep="\n")
}