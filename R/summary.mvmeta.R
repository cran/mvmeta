summary.mvmeta <-
function(object, ci.level=0.95, digits=4, ...) {

	if(ci.level<=0||ci.level>=1) stop("'ci.level' must be within 0 and 1")

	# CREATE USEFUL OBJECTS
	methodname <- c("reml","ml","fixed")
	methodlabel <- c("REML","ML","Fixed")
	p <- object$dim$p
	int <- object$int

###########################################################################
# HEADING AND SUBHEADING

	# HEADING
  cat("call:\n")
  print(object$call)
  cat("\n")

	# SUB-HEADING
  cat(if(object$dim$k==1)"Uni" else "Multi","variate ",
		ifelse(object$method=="fixed","fixed","random"),"-effects meta-",
		ifelse(p-int>0,"regression","analysis"),"\n",sep="")
	# CHECK LATER FOR CHOICE META-ANALYSIS OR METAREGRESSION
	cat("Dimension: ",object$dim$k,"\n","Studies: ",object$dim$m,"\n",sep="")
	if(object$method!="fixed") {
		cat("Estimation method: ",
			methodlabel[which(object$method==methodname)],"\n",sep="")
		cat("Variance-covariance matrix Psi: ","unstructured","\n",sep="")
	}
	cat("\n")

###########################################################################
# FIXED EFFECTS ESTIMATES

	cat("Fixed effects","\n",sep="")

	# COMPUTE USEFUL OBJECTS
	coef.se <- sqrt(diag(object$vcov))
	zval <- object$coef/coef.se
	zvalci <- qnorm((1-ci.level)/2,lower.tail=FALSE)
	pvalue <- 2*(1-pnorm(abs(zval)))
	ci.lb <- object$coef-zvalci*coef.se
	ci.ub <- object$coef+zvalci*coef.se
	cilab <- paste(round(ci.level,2)*100,"%ci.",c("lb","ub"),sep="")
  signif <- symnum(pvalue, corr = FALSE, na = FALSE, cutpoints = c(0, 
	  0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
		"*", ".", " "))

	# FOR SIMPLE META-ANALYSIS
	if(p-int==0) {
		table <- cbind(object$coef,coef.se,zval,pvalue,ci.lb,ci.ub)
		rownames(table) <- object$lab$klab
		colnames(table) <- c("Estimate","StdErr","z","p-value",cilab)
		table <- formatC(table,digits=digits,format="f")
		table <- cbind(table,signif)
		colnames(table)[7] <- ""
		print(table,quote=FALSE,right=TRUE,print.gap=2)
	}	

	# FOR META-REGRESSION
	if(p-int>0) {
		p <- object$dim$p
		tabletot <- cbind(object$coef,coef.se,zval,pvalue,ci.lb,ci.ub)
		for(i in seq(object$dim$k)) {
			ind <- seq((i-1)*p+1,(i-1)*p+p)
			table <- tabletot[ind,,drop=FALSE]
			rownames(table) <- object$lab$plab
			colnames(table) <- c("Estimate","StdErr","z","p-value",cilab)
			table <- formatC(table,digits=digits,format="f")
			table <- cbind(table,signif[ind])
			colnames(table)[7] <- ""
			cat(object$lab$klab[i],":","\n")
			print(table,quote=FALSE,right=TRUE,print.gap=2)
		}
	}
  cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")

###########################################################################
# RANDOM EFFECTS ESTIMATES

	if(!object$method=="fixed") {
	
		cat("Variance components: between-studies stdev and correlation matrix",
			"\n",sep="")

		# COMPUTE USEFUL OBJECTS
		Psisd <- sqrt(diag(object$Psi))
		Psicor <- round(cov2cor(object$Psi),3)
		Psicor[upper.tri(Psicor)] <- NA

		table <- cbind(Psisd,Psicor)
		colnames(table)[1] <- "StdDev"
		table <- formatC(table,digits=digits,format="f")
		table[grep("NA",table)] <- "."
		print(table,quote=FALSE,right=TRUE,na.print="",print.gap=2)

		cat("\n")
	}

  # OVERALL QTEST AND I-SQUARE
	q <- qtest(object)
  Q <- formatC(q$Q,digits=digits,format="f")
  pvalue <- formatC(q$pvalue,digits=digits,format="f")
  i2 <- formatC(pmax((q$Q-q$df)/q$Q*100,1),digits=1,format="f")
	cat(if(q$k==1) "Uni" else "Multi","variate ","Cochran Q-test for ",
  	if(q$residual) "residual ", "heterogeneity:","\n",sep="")
  cat("Q = ",Q[1]," (df = ",q$df[1],"), p-value = ",pvalue[1],"\n",sep="")
  cat("I-square statistic = ",i2[1],"%","\n\n",sep="")
  
###########################################################################
# FIT STATS

	cat(object$df$nobs," observations, ",object$df$fixed," fixed and ",
		object$df$random," random parameters","\n",sep="")
	table <- c(logLik(object),AIC(object),BIC(object))
	names(table) <- c("logLik","AIC","BIC")
	table <- formatC(table,digits=digits,format="f")
	print(table,quote=FALSE,right=TRUE,print.gap=2)
	cat("\n")
}

