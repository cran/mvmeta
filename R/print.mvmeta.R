print.mvmeta <-
function(x, digits=4, ...) {

	# CREATE USEFUL OBJECTS
	methodname <- c("reml","ml","fixed")
	methodlabel <- c("REML","ML","Fixed")
	p <- x$dim$p
	int <- x$int

###########################################################################
# HEADING AND SUBHEADING

  # HEADING
  cat("call:\n")
  print(x$call)
  cat("\n")

	# SUB-HEADING
  cat(if(x$dim$k==1)"Uni" else "Multi","variate ",
		ifelse(x$method=="fixed","fixed","random"),"-effects meta-",
		ifelse(p-int>0,"regression","analysis"),"\n",sep="")
	# CHECK LATER FOR CHOICE META-ANALYSIS OR METAREGRESSION
	cat("Dimension: ",x$dim$k,"\n","Studies: ",x$dim$m,"\n",sep="")
	if(x$method!="fixed") {
		cat("Estimation method: ",
			methodlabel[which(x$method==methodname)],"\n",sep="")
		cat("Variance-covariance matrix Psi: ","unstructured","\n",sep="")
	}
	cat("\n")

###########################################################################
###########################################################################
# FIT STATS

  cat(x$df$nobs," observations, ",x$df$fixed," fixed and ",
		x$df$random," random parameters","\n",sep="")
	table <- c(logLik(x),AIC(x),BIC(x))
	names(table) <- c("logLik","AIC","BIC")
	table <- formatC(table,digits=digits,format="f")
	print(table,quote=FALSE,right=TRUE,print.gap=2)
	cat("\n")

}

