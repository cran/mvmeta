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
	#cat("\n")
	cat(if(object$dim$k==1)"UNI" else "MULTI","VARIATE ",
		ifelse(object$method=="fixed","FIXED","RANDOM"),"-EFFECTS META-",
		ifelse(!is.null(object$X),"REGRESSION","ANALYSIS"),sep="")
	# CHECK LATER FOR CHOICE META-ANALYSIS OR METAREGRESSION
	cat("\n")

	# SUB-HEADING
	cat("Dimensions: ",object$dim$k,"\n","Studies: ",object$dim$m,"\n",sep="")
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
	beta.se <- sqrt(diag(object$vcov))
	zval <- object$beta/beta.se
	zvalci <- qnorm((1-ci.level)/2,lower.tail=FALSE)
	pvalue <- 2*(1-pnorm(abs(zval)))
	ci.lb <- object$beta-zvalci*beta.se
	ci.ub <- object$beta+zvalci*beta.se
	cilab <- paste(round(ci.level,2)*100,"%ci.",c("lb","ub"),sep="")
      signif <- symnum(pvalue, corr = FALSE, na = FALSE, cutpoints = c(0, 
		0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
		"*", ".", " "))

	# FOR SIMPLE META-ANALYSIS
	if(p-int==0) {
		table <- cbind(object$beta,beta.se,zval,pvalue,ci.lb,ci.ub)
		rownames(table) <- object$lab$klab
		colnames(table) <- c("Estimate","StdErr","z","p-value",cilab)
		table <- formatC(table,digits=digits,format="f")
		table <- cbind(table,signif)
		colnames(table)[7] <- ""
		print(table,quote=FALSE,right=TRUE,print.gap=2)
		cat("\n")
	}	

	# FOR META-REGRESSION
	if(p-int>0) {
		p <- object$dim$p
		tabletot <- cbind(object$beta,beta.se,zval,pvalue,ci.lb,ci.ub)
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
			cat("\n")
		}
	}

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

	q <- qtest(object)
	print(q)

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

