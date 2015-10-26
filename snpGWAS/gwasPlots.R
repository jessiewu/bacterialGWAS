## Plot manhattan plot for log reg
manhattan.plot<-function(
	position = NULL,
	pvalue = NULL,
	filename = NULL,
	plot.title = NULL){
	png(file = filename, width=1000, height=700)
	plot(x = position, y = pvalue, type="p",
		main = plot.title,
		xlab = "Position in reference genome",
		ylab = "-log10(p)"
	)
	dev.off()
}