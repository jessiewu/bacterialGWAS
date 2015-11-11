#panGWAS.R runs GWAS based on logistic regression for a pan-genome dataset where gene presence or abasence is the binary genotype.

help = paste("", "panGWAS.r runs basic and lmm GWAS using presence or absence of each gene in a pan-genome as the genotype", "Usage: Rscript panGWAS.r, -pangenome  -genomes -prefix","Options:","-pangenome","File containing matrix of gene presence/absence for each COG in a pan-genome. One row per COG, 'COGid, p/a p/a p/a'","E.g. 'gene1  1  0  0	0	0	0	1	0	0' in a tabular format","-genomes","list of genomes in sample","-phenotypes","binary phenotypes coded as 1 (case) and 0 (control) for each genome in sample","—relatedness_matrix","relatedness matrix calculated from SNPs in the dataset; output of CFML/gemma","-prefix","Output file name prefix","",sep="\n")

#read arguments: pangenome, genome list, relatedness matrix (gemma/cfml output), matrix of phenotypes (min 2 cols: genome ids + 1 column for each phenotype, with cases coded as 1 and controls coded as 0) ), output prefix
args = commandArgs(trailingOnly = TRUE)
if(length(args!=0)){
  if(args[1]=="-help" | args[1]=="-h"){
    cat(help,sep="\n")
    q("no")
  }
}

if((length(args)%%2)!=0 | length(args)==0) {
  cat(help,sep="\n")
  stop("\nNo arguments supplied\n")
}

s=seq(2,by=2,length(args))
args2=args[seq(1,by=2,length(args))]
inputs=matrix(0,ncol=2,nrow=length(args2))
for(i in 1:length(args2)){
  arg=unlist(strsplit(args2[i],""))
  inputs[i,1]=paste(arg[2:length(arg)],collapse="")
  inputs[i,2]=args[s[i]]
}
cat(args)
pg=c(0)
genomes=c(0)
ptype=c(0)
relm = c(0)
prefix=c(0)


for(i in 1:length(inputs[,1])){
  if(inputs[i,1]=="pangenome"){
    pg=inputs[i,2]
  }
  else if(inputs[i,1]=="genomes"){
    genomes=inputs[i,2]
  }
  else if(inputs[i,1]=="relatedness_matrix"){
    relm=inputs[i,2]
  }
  else if(inputs[i,1]=="phenotypes"){
    ptype=inputs[i,2]
  }
  else if(inputs[i,1]=="prefix"){
    prefix=inputs[i,2]
  }
}

existence = file.exists(pg)
if(length(which(existence==FALSE))>0) {
  stop("\n Pangenome file doesn't exist\n")
}
existence = file.exists(genomes)
if(length(which(existence==FALSE))>0) {
  stop("\n Genome list doesn't exist\n")
}
existence = file.exists(relm)
if(length(which(existence==FALSE))>0) {
  stop("\n Relatedness matrix doesn't exist\n")
}
existence = file.exists(ptype)
if(length(which(existence==FALSE))>0) {
  stop("\n Phenotypes file doesn't exist\n")
}
#load files
pangenome=read.table(pg, header=F)
genome_list = scan(genomes, what=character(0))
pheno = read.table(ptype, header=T)

#check that columns in pan-genome are ordered correctly
ids=as.character(pangenome[,1])
geno = pangenome[,2:ncol(pangenome)]
#name genomes sampled with genomes from input list
genome_list=substr(genome_list,1,9) #here I've assume 9 character comids as ids
names(geno)=genome_list
#order genomes in pan-genome matrix to match order of genomes in relatedness matrix/phenotype file (/dipro/mmm/gorm/v3/ana/Saur/AbxGWAS/derval_gwas/der_and_val_paths_phenotype.txt)
#specify to user that genome ids need to be the same (but not in the same order) for the list of genomes input to create the pan-genome and for the genome ids used in the phenotype file (ie only file suffixes differ)
geno = geno[,match(pheno$Comid, names(geno))]
pg_ordered = cbind(ids,geno)
write.table(pg_ordered, paste0(pg,".ordered"), row.names=F,col.names=F,quote=F,sep="\t")

#create gemma_gen_format.txt and gemma_snp_format.txt files for lmm GWAS
#gemma_gen_format = ids, 1, 0, genotypes (1/0)
gen=cbind(ids, rep(1,nrow(pangenome)),rep(0,nrow(pangenome)),geno)
write.table(gen,paste0(prefix,".gemma_gen.txt"),row.names=F, col.names=F, quote=F)
#gemma_snp_format = snp_id, position, chr (24 for haploids)
snp=cbind(ids,seq(1:nrow(pangenome)),rep(24, nrow(pangenome)))
write.table(snp,paste0(prefix,".gemma_snp.txt"),row.names=F, col.names=F, quote=F)

#run vanilla and lmm GWAS for given set of phenotypes and pan-genome set
getLogNull=function(datafiles){ #function for correcting p-values output from gemma
  a=scan(datafiles,what=character(0),sep="\n")
  b=unlist(strsplit(a[17]," "))
  b=as.numeric(b[length(b)])
  return(b)
}
for (p in 2:ncol(pheno)){ #in this file the first 3 columns are additional info: /dipro/mmm/gorm/v3/ana/Saur/AbxGWAS/derval_gwas/der_and_val_paths_phenotype.txt
  #create phenotype files
  name=names(pheno)[p]
  cat(paste0("Doing: ",name,"\n"))
  write.table(pheno[,p],paste0(name,".pheno"), row.names=F, col.names=F,quote=F)
  #vanilla GWAS
  pt = paste0(name,".pheno")
  outfix = paste0(name,"_",prefix)
  gwas = paste0("Rscript /dipro/mmm/gorm/v2/ana/Saur/jchar/pangenome/code/panGWAS/do_logreg_chw.R -bips ",paste0(pg,".ordered")," -phenotype ",pt," -prefix ",paste0(outfix,".raw"))
  cat(paste0(gwas,"\n"))
  system(gwas)
  #lmm GWAS
  #run gemma command *this path works*
  lmm_gwas=paste0("/dipro/mmm/gorm/v3/ana/Cdif/earle/SNPanalysis/gemma/gemma -g ",paste0(prefix,".gemma_gen.txt")," -p ",pt," -a ",paste0(prefix,".gemma_snp.txt")," -k ",relm," -lmm 4 "," -o ",paste0(outfix,".lmm")," -maf 0")
  cat(paste0(lmm_gwas,"\n"))
  system(lmm_gwas)
  #correct p-values and calculate -logl0 pvalues ###needs checking###
  logNull=getLogNull(paste0("./output/",outfix,".lmm.log.txt"))
  gem=read.table(paste0("./output/",outfix,".lmm.assoc.txt"),header=TRUE,as.is=TRUE)
  getPvalue=function(gem,logNull){
    D=2*abs(logNull-as.numeric(gem[length(gem)])) # This is the LH1 column
    pval=-log10(exp(1))*pchisq(D,1,low=F,log=TRUE)
    return(pval)
  }
  negLog10=apply(gem,1,getPvalue,logNull=logNull)
  gem_all=cbind(gem,negLog10)
  write.table(gem_all,file=paste0(outfix,"_LMM_corrected.txt"),row=F,col=T,sep="\t",quote=FALSE)
}
  
