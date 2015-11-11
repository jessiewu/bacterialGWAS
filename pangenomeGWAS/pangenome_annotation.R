# Annotates genes in pan-genome using BLAST against both a refseq database plus the nucleotide database
# Jane Charlesworth 2014

## pan-genome annotation script
help=paste("","pangenome_annotation.R","Usage: Rscript pangenome_annotation.R -pangenome -all_genes -prefix -blastdb1 -blastdb2 -ncbi_db -nproc","       E.g. Rscript pangenome_annotation.R -pangenome prefix.cdhit.varsites -all_genes all_prot.faa -present_ctrl -prefix annotationrun1 -blastdb1 db1 -blastdb2 db2 -ncbi_db genedb.txt -nproc 8 ","","Options:","-pangenome"," Matrix of 1s and 0s denoting presence or absence of each variable gene in pan-genome","-all_genes"," FASTA file containing all annotated coding sequences used to construct pan-genome","-prefix"," File name prefix for output files","-blastdb1"," First blast database to use to annotate kmers","-blastdb2"," Second blast database to use to annotate kmers (more extensive search)","-ncbi_db"," Gene database to annotate blast results","-nproc"," Number of proccessors to run blast in parallel (Default = 1)",sep="\n")

library(seqinr)

#######################################################################################

## Functions

#######################################################################################

#######################################################################################
## Read in genes in pan-genome
## Opens file of all coding sequences ("all_prot.faa") and pan-genome file from a CD-hit
## run and creates a fasta file containing only the protein sequences in a given
## pan-genome created by CD-hit.
#######################################################################################
read_genes=function(all_genes,pangenome){
  ## Read in set of all CDS
  cds = read.fasta(all_genes, seqtype="AA")
  names = attr(cds,"name")
  ## Read in pan-genome
  pg=read.table(pangenome, header=T, stringsAsFactors=F)
  pg_names=pg[,1]
  pg_names=substr(pg_names,2,nchar(pg_names)) #trim > character from names
  seqs=cds[which(pg_names %in% names)]
  return(seqs)
}
#######################################################################################
## Split the genes into chunks of 1000 for faster annotation
##Unlike with the kmers, all genes in pan-genome need to be annotated, so don't worry
##about their significance
#######################################################################################
split_genes=function(seqs){
	nsplit=floor(length(seqs)/1000)
	splitseq=seq(from=1,by=1000,to=1000*nsplit)
	splitseq2=c(splitseq[2:length(splitseq)]-1,length(seqs))
	return(cbind(splitseq,splitseq2))
}	

#######################################################################################
## Construct a file containing the genes for the BLAST search
##From this point onwards, gene1000 is the list of names of the 1000 genes in each set
## of genes input to BLAST
#######################################################################################
construct_gene_file=function(data_loop,seqs){
  gene1000=attr(seqs[data_loop],"name")
  write.fasta(seqs[gene1000],gene1000,paste0(prefix,"temp_genes.txt"))
}

#######################################################################################
## Command to run BLAST
#######################################################################################	
run_BLAST=function(prefix,nproc,db){
	system(paste0("tblastn -query ",prefix,"temp_genes.txt -db ",db," -evalue 10 -num_threads ",nproc," -outfmt '6 qseqid sseqid sacc pident length mismatch gapopen qstart qend evalue sstart send qseq' -max_target_seqs 1 -out ",prefix,"_blast.txt"))
}


#######################################################################################
## Read in the BLAST result	
#######################################################################################	
read_blast_result=function(blastfile){
	blast.search=scan(blastfile,what=character(0))
	blast.search=matrix(blast.search,ncol=13,byrow=TRUE)
}

#######################################################################################
## Check which genes had a result in the blast search
## check needs names of genes from each set of 1000
#######################################################################################
gene_check=function(gene1000,blast.search){
	check=length(which(blast.search[,1]==gene1000))
	return(check)
}

#######################################################################################
## If more than the top hit was returned, remove any lower hits
#######################################################################################
remove_multiple=function(which.multiples,blast.search,gene1000){
	blast.multiple=which(blast.search[,1]==gene1000[which.multiples])
	toremove=blast.multiple[2:length(blast.multiple)]
}


#######################################################################################
## Checks which genes had a result and which rows need to be removed
## Uses functions gene_check and remove_multiple
#######################################################################################
check_genes_remove_duplicates=function(gene1000,blast.search){
	gene1000check=sapply(gene1000,gene_check,blast.search=blast.search,USE.NAMES=FALSE)
	which.multiples=which(gene1000check>1)
	toremove=unlist(sapply(which.multiples,remove_multiple,blast.search=blast.search,gene1000=attr(gene1000,"name")))
	if(length(toremove)!=0){
		blast.search=blast.search[-toremove,]
	}
	return(blast.search)
}

#######################################################################################
## Get the accession numbers and positions from the BLAST results	
#######################################################################################	
get_blast_accessions=function(blast.search){
	# Accession numbers and positions
	blast.gen=blast.search[2]
	# E-value
	blast.e=as.numeric(blast.search[10])
	# Start position
	blast.r1=as.numeric(blast.search[11])
	# End position
	blast.r2=as.numeric(blast.search[12])
	# Percentage cover
	percover=as.numeric(blast.search[5])
	return(c(blast.gen,blast.e,blast.r1,blast.r2,percover))
}	

#######################################################################################
## Function to find entry with NC_ and position from BLAST
#######################################################################################

annotation_index = function(pos,ncbi.pos){
	sapply(pos, function(i) which(ncbi.pos[,1] <=i & ncbi.pos[,2]>=i))
}
annotation_index2 = function(pos,ncbi.pos){
	sapply(pos, function(i) which(ncbi.pos[,2] <=i & ncbi.pos[,1]>=i))
}

#######################################################################################
## Get the shortened accession number from the longer entry
#######################################################################################	
get_short_accession=function(blast.gen){
	blast.gen.split=unlist(strsplit(blast.gen,""))
	gen.w=which(blast.gen.split=="|")
	nc=paste(blast.gen.split[(gen.w[3]+1):(gen.w[4]-1)],collapse="")
}

#######################################################################################
## For a gene matching to the gene database, gets the positions of the gene match from
##  the end of the accession number
## blast.gen: long accession entry from BLAST result, containing start and end
##  positions of the gene entry
#######################################################################################

get_positions=function(blast.gen){
	
	blast.gen.split=unlist(strsplit(blast.gen,""))
	gen.w=which(blast.gen.split=="|")
	nc=paste(blast.gen.split[(gen.w[3]+1):(gen.w[4]-1)],collapse="")
	
	# Does the accession entry include a c, if it does check where it is
	includes_c=length(which(blast.gen.split[(gen.w[4]+2):length(blast.gen.split)]=="c"))
	plus_2_c=blast.gen.split[(gen.w[4]+2)]!="c"
	# Are there only two positions in this entry
	two_positions=length(which(blast.gen.split[gen.w[4]:length(blast.gen.split)]=="-"))==1
	if(is.na(includes_c)==FALSE&is.na(plus_2_c)==FALSE&length(which(blast.gen.split[gen.w[4]:length(blast.gen.split)]=="-"))!=0){
	
	if(includes_c<=1 & plus_2_c==TRUE & two_positions==TRUE){
		pos=as.numeric(paste(blast.gen.split[(gen.w[4]+2):(which(blast.gen.split=="-")-1)],collapse=""))
		pos2=as.numeric(paste(blast.gen.split[(which(blast.gen.split=="-")[1]+1):length(blast.gen.split)],collapse=""))
	} else if(includes_c<=1 & plus_2_c==FALSE & two_positions==TRUE){
		pos=as.numeric(paste(blast.gen.split[(gen.w[4]+3):(which(blast.gen.split=="-")-1)],collapse=""))
		pos2=as.numeric(paste(blast.gen.split[((which(blast.gen.split=="-")+1)):length(blast.gen.split)],collapse=""))
	} else if(includes_c<=1 & plus_2_c==FALSE & two_positions==FALSE){
		pos1=as.numeric(paste(blast.gen.split[(gen.w[4]+3):(which(blast.gen.split=="-")[1]-1)],collapse=""))
		pos2=as.numeric(paste(blast.gen.split[((which(blast.gen.split=="-")[1]+1)):(gen.w[4]+(which(blast.gen.split[(gen.w[4]+2):length(blast.gen.split)]=="c")[2])-1)],collapse=""))
		pos3=as.numeric(paste(blast.gen.split[(gen.w[4]+which(blast.gen.split[(gen.w[4]+2):length(blast.gen.split)]=="c")[2]+2):(which(blast.gen.split=="-")[2]-1)],collapse=""))
		pos4=as.numeric(paste(blast.gen.split[(which(blast.gen.split=="-")[2]+1):length(blast.gen.split)],collapse=""))
				
		pos=min(pos1,pos2,pos3,pos4)
		pos2=max(pos1,pos2,pos3,pos4)
		rm(pos1)
		rm(pos3)
		rm(pos4)
		
	} else if(includes_c==0 & two_positions==FALSE){
		pos1=as.numeric(paste(blast.gen.split[(gen.w[4]+2):(which(blast.gen.split=="-")[1]-1)],collapse=""))
		pos2=as.numeric(paste(blast.gen.split[((which(blast.gen.split=="-")[1]+1)):(gen.w[4]+which(blast.gen.split[gen.w[4]:length(blast.gen.split)]==",")-2)],collapse=""))
		pos3=as.numeric(paste(blast.gen.split[(gen.w[4]+which(blast.gen.split[gen.w[4]:length(blast.gen.split)]==",")):(which(blast.gen.split=="-")[2]-1)],collapse=""))
		pos4=as.numeric(paste(blast.gen.split[(which(blast.gen.split=="-")[2]+1):length(blast.gen.split)],collapse=""))
				
		pos=min(pos1,pos2,pos3,pos4)
		pos2=max(pos1,pos2,pos3,pos4)
		rm(pos1)
		rm(pos3)
		rm(pos4)
	} else {
		pos=0
		pos2=0
	}
	
} else {
	pos=0
	pos2=0
}
return(c(pos,pos2,nc))
}


#######################################################################################
## Get the annotations from the gene database
## Uses the function annotation_index
## positions: a 4 column matrix of the start and end positions of each gene, accession
##  number and beginning position of the gene
## ncbi.gene: the gene database
#######################################################################################

get_annotations=function(positions,ncbi.gene){
	
	pos=as.numeric(positions[1])
	pos2=as.numeric(positions[2])
	nc=positions[3]
	blast.r1=as.numeric(positions[4])
	
	# If there are matches to the accession number in the gene database
	# Get the rows of the gene database which match the accession
	if(is.na(pos)==FALSE & is.na(pos2)==FALSE & length(which(ncbi.gene[,6]==nc))!=0){
		ncbi.match=ncbi.gene[which(ncbi.gene[,6]==nc),]
		# Get the positions of the genes of these matches
		# If ncbi.match gives just one match, check 
		if(length(ncbi.match)==9){
			ncbi.pos=matrix(as.numeric(ncbi.match[7:8]),ncol=2)
		} else {
			ncbi.pos=matrix(as.numeric(unlist(ncbi.match[,7:8])),ncol=2)
		}
			
		# Find which gene out of the selection with matching accession number based on position
		if(pos<pos2){
			annot=unlist(sapply((pos+blast.r1),annotation_index,ncbi.pos=ncbi.pos))[1]
		} else {
			annot=unlist(sapply((pos2+blast.r1),annotation_index,ncbi.pos=ncbi.pos))[1]
		}

		# If the positions match for one of the entries
		if(length(annot)!=0){
			if(annot!="integer(0)" & is.na(annot)==FALSE){
				if(length(ncbi.match)==9){
					annot2=ncbi.match
				} else {
					annot2=ncbi.match[annot,]
				}
				
			} else if(annot=="integer(0)" & is.na(annot)==FALSE){
				annot=annotation2(pos+blast.r1)[1]
				if(length(ncbi.match)==9){
					annot2=ncbi.match
				} else {
					annot2=ncbi.match[annot,]
				}
			} else if(is.na(annot)==TRUE){
				annot2=matrix("NA",9,nrow=1)
			}
		} else {
			annot2=matrix("NA",9,nrow=1)
		}
	} else {
		annot2=matrix("NA",9,nrow=1)
	}		
}


#######################################################################################
## Function to call the other functions to annotate the genes
##throughout this "gene1000" = list of names of genes in each BLAST input set
#######################################################################################	
annotate_genes=function(splitseq,prefix,nproc,blastdb1,blastdb2,seqs,gene_output_file){
	data_loop=seq(from=splitseq[1],by=1,to=splitseq[2])
	
	### Annotate using the genus BLAST database and specified gene database
	# Construct file containing 1000 genes
	construct_gene_file(data_loop,seqs=seqs)
	gene1000=attr(seqs[data_loop],"name") #make list of names for 1000 genes in each set
	## Run BLAST from the command line against refseq database and store wanted values
	run_BLAST(prefix,nproc,blastdb1)
	# Open BLAST results
	blastfile=paste0(prefix,"_blast.txt")
	blast.search=read_blast_result(blastfile)
	# gene check - find which of the 1000 genes have a result
	blast.search=check_genes_remove_duplicates(gene1000,blast.search)
	# Which genes had no BLAST result		
	check_noresult=which(sapply(gene1000,gene_check,blast.search=blast.search,USE.NAMES=FALSE)==0)
	# Which genes did have a BLAST result	
	whichrun=setdiff(1:1000,check_noresult)
	index=whichrun
	# Get BLAST accessions
	blast_acc=apply(blast.search,1,get_blast_accessions)
	blast_acc=matrix(blast_acc,ncol=5,byrow=TRUE)
	# Get positions
	positions=sapply(blast_acc[,1],get_positions)
	positions=matrix(positions,ncol=3,byrow=TRUE)
	# Get annotations for those with a BLAST result to the genus BLAST database
	annotation=apply(cbind(positions,blast_acc[,3]),1,get_annotations,ncbi.gene=ncbi.gene)
	annotation=matrix(annotation,ncol=9,byrow=TRUE)
	info1=cbind(whichrun,gene1000[data_loop[1]+whichrun-1],blast_acc,positions[,3],annotation)
	# When there is no match to the genus BLAST database
	if(length(check_noresult)>0){
		annotation=matrix("NA",15,nrow=length(check_noresult))
		info2=cbind(check_noresult,gene1000[data_loop[1]+check_noresult-1],annotation)
		}
	if(length(check_noresult)!=0){
		info.a=rbind(info1,info2)
	} else {
		info.a=info1
	}
	
	#### Also run BLAST on second database, the whole nucleotide database
	# Run BLAST	
	run_BLAST(prefix,nproc,blastdb2)
	# Open BLAST results
	blastfile=paste0(prefix,"_blast.txt")
	blast.search=read_blast_result(blastfile)
	# gene check - find which of the 1000 genes have a result
	blast.search=check_genes_remove_duplicates(data_loop,blast.search)
	# Which genes had no BLAST result		
	check_noresult=which(sapply(gene1000,gene_check,blast.search=blast.search,USE.NAMES=FALSE)==0)
	# Which genes did have a BLAST result	
	whichrun=setdiff(1:1000,check_noresult)
	index=whichrun
	# Get BLAST accessions
	blast_acc=apply(blast.search,1,get_blast_accessions)
	blast_acc=matrix(blast_acc,ncol=5,byrow=TRUE)
	# Get the short accession number
	nc=sapply(blast_acc[,1],get_short_accession,USE.NAMES=FALSE)
	info3=cbind(whichrun,gene1000[data_loop[1]+whichrun-1],blast_acc,nc)	
	if(length(check_noresult)>0){
		annotation=matrix("NA",6,nrow=length(check_noresult))
		info4=cbind(check_noresult,gene1000[data_loop[1]+check_noresult-1],annotation)
	}
	
	if(length(check_noresult)!=0){
		info.b=rbind(info3,info4)
	} else {
		info.b=info3
	}
	info.all=cbind(info.a,info.b)

	write.table(info.all,file=gene_output_file,quote=F,row.names=F,col.names=F,sep="\t",append=TRUE)
	cat(paste0("Annotated ",data_loop[length(data_loop)]," genes ",Sys.time()),file=paste0(prefix,"_annotlogfile.txt"),sep="\n",append=TRUE)
	cat("Done",data_loop[length(data_loop)],"\n")
	
}




#######################################################################################
## Read in command line arguments
#######################################################################################

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

# Separate command line arguments from their -names
s=seq(2,by=2,length(args))
args2=args[seq(1,by=2,length(args))]
inputs=matrix(0,ncol=2,nrow=length(args2))
for(i in 1:length(args2)){
	arg=unlist(strsplit(args2[i],""))
	inputs[i,1]=paste(arg[2:length(arg)],collapse="")
	inputs[i,2]=args[s[i]]
}

pangenome=c()
all_genes=c()
prefix=c()
blastdb1=c()
blastdb2=c()
ncbi_db=c()
nproc=as.numeric(1)

for(i in 1:length(inputs[,1])){
  if(inputs[i,1]=="prefix"){
		prefix=inputs[i,2]
		}
  else if(inputs[i,1]=="all_genes"){
    all_genes = inputs[i,2]
  }
  else if(inputs[i,1]=="pangenome"){
    pangenome = inputs[i,2]
  }
	else if(inputs[i,1]=="blastdb1"){
		blastdb1 = inputs[i,2]
		}
	else if(inputs[i,1]=="blastdb2"){
		blastdb2 = inputs[i,2]
		}
	else if(inputs[i,1]=="ncbi_db"){
		ncbi_db=inputs[i,2]
	}
	else if(inputs[i,1]=="nproc"){
		nproc=as.numeric(inputs[i,2])
	}
}

#######################################################################################
## Check files exist and basic input checks
#######################################################################################

check_existence=c(all_genes,pangenome,blastdb1,blastdb2,ncbi_db)


for(i in 1:length(check_existence)){
	check=check_existence[i]
	if(check!=blastdb1 & check!=blastdb2){
		if(file.exists(check)==FALSE){
			stop("\nIncorrect usage: ",check," file doesn't exist\n")
		}
	} else {
		dbpath=unlist(strsplit(check,"/"))
		db=dbpath[length(dbpath)]
		dbpath=paste(dbpath[1:(length(dbpath)-1)],collapse="/")
		existence=dir(dbpath,pattern=glob2rx(paste0(db,"*")))
		if(length(existence)<3){
			stop("\nIncorrect usage: ",check," database doesn't exist")
		}
	}
}


if(nproc>as.numeric(system("nproc",intern=TRUE))){
	stop("\nIncorrect usage: Number of processors specified is greater than the number of the machine")
}
if((nproc%%1)!=0 | nproc==0){
	stop("\nIncorrect usage: number of processors must be an integer\n")
}



#######################################################################################
## Read in files
######################################################################################
#read in gene file and get just genes in pan-genome
seqs= read_genes(all_genes,pangenome)

# Split into chunks of 1000 to BLAST simultaneously - much faster
splitseq=split_genes(seqs)


# Read in gene database
ncbi.gene=scan(paste0(ncbi_db),what=character(0),sep="\t",quote="")
ncbi.gene=matrix(ncbi.gene,ncol=9,byrow=TRUE)


#######################################################################################
# Output file name and header information
#######################################################################################

## Output file name
gene_output_file=paste0(prefix,".out.annotated.txt")
# Write the annotation start time to a log file
cat(paste0("Annotation start time ",Sys.time()),file=paste0(prefix,"_annotlogfile.txt"),sep="\n",append=TRUE)
# Write headers to file
#headerinfo=c("Gene","Accession_all","Evalue","Start_pos","End_pos","Length_match","Accession","Gene_db_index","Gene","Product","Other_aliases","Genomic_context","Accession","Gene_start","Gene_end","Gene_starta","Gene_endb","Gene_startc","Gene_endc","Other_designations","Gene","Accession_all","Evalue","Start_pos","End_pos","Length_match","Accession")
headerinfo=c("Gene","Accession_all","Evalue","Start_pos","End_pos","Length_match","Accession","Gene_db_index","Gene","Product","Other_aliases","Genomic_context","Accession","Gene_start","Gene_end","Other_designations","Gene","Accession_all","Evalue","Start_pos","End_pos","Length_match","Accession")
cat(paste(c(headerinfo),collapse="\t"),file=gene_output_file,append=FALSE,sep="\n")

#######################################################################################
## Run annotation calling annotate_gene function
#######################################################################################

apply(splitseq,1,annotate_genes,prefix=prefix,nproc=nproc,blastdb1=blastdb1,blastdb2=blastdb2,seqs=seqs,gene_output_file=gene_output_file)

cat(paste0("Annotation of genes ended: ",Sys.time()),file=paste0(prefix,"_annotlogfile.txt"),append=TRUE,sep="\n")
