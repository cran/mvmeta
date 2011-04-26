print.qtest.mvmeta <-
function(x, digits=4, ...) {

	if(digits<1) stop("'digits' must be >=1")
	Q <- formatC(x$Q,digits=digits,format="f")
	if(x$pvalue<10^-digits) {
		pvalue <- paste("<","0.",paste(rep("0",digits-1),collapse=""),
			1,sep="")
	} else pvalue <- formatC(x$pvalue,digits=digits,format="f")

	cat(if(x$k==1) "Uni" else "Multi","variate ","Cochran Q-test for ",
		if(x$residual) "residual ", "heterogeneity:","\n",sep="")
	cat("Q = ",Q," (df = ",x$df,"), p-value = ",pvalue,"\n\n",sep="")
}

