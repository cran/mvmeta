print.qtest.mvmeta <-
function(x, digits=3, ...) {

	if(digits<1) stop("'digits' must be >=1")
  signif <- symnum(x$pvalue, corr = FALSE, na = FALSE, cutpoints = c(0, 
  	0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
		"*", ".", " "))
  Q <- formatC(x$Q,digits=digits,format="f")
  pvalue <- formatC(x$pvalue,digits=digits,format="f")
  table <- cbind(Q,df=x$df,"p-value"=pvalue,signif)
  colnames(table)[4] <- ""

	cat(if(x$k==1) "Uni" else "Multi","variate ","Cochran Q-test for ",
		if(x$residual) "residual ", "heterogeneity","\n",sep="")
  if(x$k>1) cat("\nOverall test:","\n")
	cat("Q = ",Q[1]," (df = ",x$df[1],"), p-value = ",pvalue[1],"\n",sep="")
  if(x$k>1) cat("\nTests on single outcome parameters:","\n")
  if(x$k>1) print(table[-1,,drop=FALSE],quote=FALSE,right=TRUE,print.gap=2)
  if(x$k>1) cat("---\nSignif. codes: ", attr(signif, "legend"), "\n")
  cat("\n")
}

