
species = c("Saur","Ecol","Kleb","Mtub")
par = c(c(0.7,0.7),c(0.8,0.8),c(0.5,0.5))

library(seqinr)
#read arguments: cdhit input, species, identity, coverage, wordsize, whether to run cdhit or not (t/f),whether to run cdhit or not (t/f), list of annotated genomes
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
input=c(0)
species=c(0)
id=c(0)
cov=c(0)
wsize=c(0)
cdhit=c(0)
confile=c(0)


for(i in 1:length(inputs[,1])){
  if(inputs[i,1]=="input"){
    input=inputs[i,2]
  }
  else if(inputs[i,1]=="id"){
    id=inputs[i,2]
  }
  else if(inputs[i,1]=="cov"){
    cov=inputs[i,2]
  }
  else if(inputs[i,1]=="species"){
    species=inputs[i,2]
  }
  else if(inputs[i,1]=="wsize"){
    wsize=inputs[i,2]
  }
  else if(inputs[i,1]=="cdhit"){
    cdhit=inputs[i,2]
  }
  else if(inputs[i,1]=="contig_files"){
    confile=inputs[i,2]
  }
}
existence = file.exists(input)
if(length(which(existence==FALSE))>0) {
  stop("\n Input file doesn't exist\n")
}

#######################################
outfix = paste0(species,".",id,".",cov)
if (cdhit=="t"){ #code to rename contigs within assembly files (eg ecol and kleb - needs testing)
  system(paste0("cd-hit -i ",input," -o ",outfix," -c ",id," -M 0 -G 1 -aL ",cov," -aS ",cov," -n ",wsize," -g 1"))
}

#read and process clusters
clusterfile = paste0(species,".",id,".",cov,".clstr")
cn <- paste("V",1:5,sep="")
clusters = read.table(clusterfile, fill=T, stringsAsFactors=F,col.names=cn)
is.newclust=clusters[,1]==">Cluster"
#get names of cluster representatives to use as ids for pan-genome
clusterrep=which(clusters[,4]=="*")
clusternames=substr(clusters[clusterrep,]$V3,1,nchar(clusters[clusterrep,]$V3)-3)
clusterid = cumsum(is.newclust)
clusterid = clusterid[clusters[,1]!=">Cluster"]
clusterid=as.factor(clusterid)
genomeid=substr(clusters[,3],2,10)
genomeid=genomeid[clusters[,1]!=">Cluster"]
tb = table(genomeid,clusterid)
genomeid=as.factor(genomeid)
write.table(tb,paste0(outfix,".cdhit.pa"), row.names=F,col.names=F,quote=F,sep="\t")
present = colSums(tb > 0)
png(paste0(outfix,".geneFreqs.png"))
#hist(present,100)
plot.ecdf(present)
title=(main=paste0(outfix,".geneFreqs"))
dev.off()

######################################
#prepare second matrix for input to GWAS code; this one is a matrix of 1s and 0s with an id column and only includes variable genes. Rows are genes in pan-genome and columns are genomes.
pa = table(clusterid,genomeid)
row.names(pa)=clusternames
pa[pa>1]=1 #recode any clusters with multiple hits in a given genome to just count presence once
vs = pa[which(rowSums(pa)!=ncol(pa)),]#get variable sites
write.table(vs,paste0(outfix,".cdhit.varSites"), row.names=T,col.names=F,quote=F,sep="\t")
write.table(levels(genomeid), paste0(outfix,".genomes"),row.names=F,col.names=F,quote=F,sep="\t") #output list of genomes for GWAS sorting


######################################
#Danny's c code to look at distribution of homologs between and within genomes
dyn.load("/home/wilson/C++/R/forJane/forJane.so")
# X[i,j] is the number of ORFs in genome i, cluster j
X.row1 = scan(paste0(outfix,".cdhit.pa"),nlines=1)
X = scan(paste0(outfix,".cdhit.pa"))
X = matrix(X,byrow=T,ncol=length(X.row1))

# Maximum number of homologs
kmax = max(X) - 1
# Number of ORFs per genome
L = rowSums(X)
# Frequency distribution for the number of homologs in the same genome, taken over all ORFs and genomes
# Y[k] is the frequency of ORFs that have k homologs in the same genome
system.time((Y = sapply(0:(kmax),function(k) mean(rowSums((X==k+1)*X)/L))))
Y = ncol(X)*Y
Y1 = vector()
cats= seq(1,length(Y))
for (i in 1:length(Y)){
  Y1 = c(Y1, rep(cats[i],Y[i]))
}
xlab=c(1:10)
png(paste0(outfix,".HomologsInGenome.png"))
barplot(Y1)
title(main=paste0(outfix,".HomologsInGenome.png"))
dev.off()
# Frequency distribution for the number of homologs in a different genome, taken over all ORFS and genomes
# Z[k] is the frequency of ORFs that have k homologs in a different genome
calcZ = function(X) {
  n = as.integer(nrow(X))
  N = as.integer(ncol(X))
  L = as.double(rowSums(X))
  kmax = as.integer(max(X))
  
  .C("calcZ",n,N,as.double(t(X)),L,kmax,ret=double(kmax+1),dup=FALSE)$ret
}

system.time((Z = calcZ(X)))
Z = ncol(X)*Z
Z1 = vector()
cats= seq(0,length(Z)-1)
for (i in 1:length(Z)){
  Z1 = c(Z1,rep(cats[i],Z[i]))
}
xlab=c(0:10)
png(paste0(outfix,".SharedHomologs.png"))
barplot(Z1)
title(main=paste0(outfix,".SharedHomologs.png"))
dev.off()

