localMinima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  
  return(y)
}

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}


getCutoff = function(filePath = NULL, prefix = NULL, width = 5){
	#message(filePath)
	counts = scan(filePath, what=integer(0),sep=",")

	if(counts[1] != 0){
		stop("The first count should be zero.")
	}
	
	counts = counts[-1]
	tab = cbind(1:length(counts),counts)
	x = unlist(apply(tab,1,function(z){rep(z[1],z[2])}))
	nbins = ceiling(length(counts)/width)
	den = density(x, n = nbins)

	
	lmIndex = localMinima(den$y)
	cutoff = floor(den$x[lmIndex ][which(den$x[lmIndex] >= 5)[1]])
	sumCounts = sum(counts)
	localMax = max(counts[cutoff:length(counts)])/sumCounts
	binSize = length(counts)/nbins
	estFreq = den$y/sum(den$y)/binSize
	
	png(file=paste(prefix,"density.png", sep="_"), width=600, height=600, res=150)
	par(mar=c(5,4,1,2)+0.2)
	plot(1:length(counts), counts/sumCounts, type="l", ylim=c(0, localMax),
	xlab="Counts", ylab="Frequency")
	lines(den$x, estFreq, col="lightblue", lty=2)
	abline(v=cutoff, col='red')
	dev.off()
	
	png(file=paste(prefix,"density2.png", sep="_"), width=600, height=600, res=150)
	par(mar=c(5,4,1,2)+0.2)
	plot(1:length(counts), counts/sumCounts, type="l",
	xlim=c(0,300), ylim=c(0, localMax),
	xlab="Counts", ylab="Frequency")
	lines(den$x, estFreq, col="lightblue", lty=2)
	abline(v=cutoff, col='red')
	dev.off()
	
	
	png(file=paste(prefix,"density3.png", sep="_"), width=600, height=600, res=150)
	par(mar=c(5,4,1,2)+0.2)
	plot(1:length(counts), counts/sumCounts, type="l",
	xlim=c(0,100), ylim=c(0, localMax),
	xlab="Counts", ylab="Frequency")
	lines(den$x, estFreq, col="lightblue", lty=2)
	abline(v=cutoff, col='red')
	dev.off()

	
	return(cutoff)
}



args = commandArgs(trailingOnly = TRUE)
if(!(length(args) == 2 | length(args) == 3| length(args) ==0)) {
	stop("\nIncorrect usage\n")
}

input.df = read.table(file=args[1], header=T, as.is=T)
prefix = args[2]
fileCount = nrow(input.df)
cutoffs = vector(length=fileCount)

width = 5
if(length(args) == 3){
	width = as.numeric(args[3])
}
message("bin width: ", width)
for(i in 1:fileCount){
	cutoffs[i] = getCutoff(filePath = input.df$filePath[i], prefix = input.df$id[i], width = width)
	cutoffs[i] = max(cutoffs[i], 5)
	
	if((i%%100) == 0){
		message(i, "complete.")
	}
}

output.df = data.frame(input.df,cutoffs)
outputFile = paste(prefix,"KmerCountCutoff.txt", sep="")
write.table(output.df, file=outputFile, row.names=F, col.names=T, quote=F, sep="\t")




