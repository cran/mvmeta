print.mvmeta <-
function(x, ci.level=0.95, digits=4, ...) {

	if(ci.level<=0||ci.level>=1) stop("'ci.level' must be within 0 and 1")

	# CREATE USEFUL xS
	methodname <- c("reml","ml","fixed")
	methodlabel <- c("REML","ML","Fixed")

###########################################################################
# HEADING AND SUBHEADING

	# HEADING
	#cat("\n")
	cat(if(x$dim$k==1)"UNI" else "MULTI","VARIATE ",
		ifelse(x$method=="fixed","FIXED","RANDOM"),"-EFFECTS META-",
		ifelse(!is.null(x$X),"REGRESSION","ANALYSIS"),sep="")
	# CHECK LATER FOR CHOICE META-ANALYSIS OR METAREGRESSION
	cat("\n\n")

###########################################################################
# FIXED EFFECTS ESTIMATES

	cat("Fixed effects","\n",sep="")

	# COMPUTE USEFUL xS
	beta.se <- sqrt(diag(x$vcov))
	zval <- x$beta/beta.se
	zvalci <- qnorm((1-ci.level)/2,lower.tail=FALSE)
	pvalue <- 2*(1-pnorm(abs(zval)))
	ci.lb <- x$beta-zvalci*beta.se
	ci.ub <- x$beta+zvalci*beta.se
	cilab <- paste(round(ci.level,2)*100,"%ci.",c("lb","ub"),sep="")
      signif <- symnum(pvalue, corr = FALSE, na = FALSE, cutpoints = c(0, 
		0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
		"*", ".", " "))

	# FOR SIMPLE META-ANALYSIS
	if(is.null(x$X)) {
		table <- cbind(x$beta,beta.se,zval,pvalue,ci.lb,ci.ub)
		rownames(table) <- x$lab$klab
		colnames(table) <- c("Estimate","StdErr","z","p-value",cilab)
		table <- formatC(table,digits=digits,format="f")
		table <- cbind(table,signif)
		colnames(table)[7] <- ""
		print(table,quote=FALSE,right=TRUE,print.gap=2)
		cat("\n")
	}	

	# FOR META-REGRESSION
	if(!is.null(x$X)) {
		p <- x$dim$p
		tabletot <- cbind(x$beta,beta.se,zval,pvalue,ci.lb,ci.ub)
		for(i in seq(x$dim$k)) {
			ind <- seq((i-1)*p+1,(i-1)*p+p)
			table <- tabletot[ind,,drop=FALSE]
			rownames(table) <- x$lab$plab
			colnames(table) <- c("Estimate","StdErr","z","p-value",cilab)
			table <- formatC(table,digits=digits,format="f")
			table <- cbind(table,signif[ind])
			colnames(table)[7] <- ""
			cat(x$lab$klab[i],":","\n")
			print(table,quote=FALSE,right=TRUE,print.gap=2)
			cat("\n")
		}
	}

###########################################################################
# RANDOM EFFECTS ESTIMATES

	if(!x$method=="fixed") {
	
		cat("Variance components: between-studies stdev and correlation matrix",
			"\n",sep="")

		# COMPUTE USEFUL xS
		Psisd <- sqrt(diag(x$Psi))
		Psicor <- round(cov2cor(x$Psi),3)
		Psicor[upper.tri(Psicor)] <- NA

		table <- cbind(Psisd,Psicor)
		colnames(table)[1] <- "StdDev"
		table <- formatC(table,digits=digits,format="f")
		table[grep("NA",table)] <- "."
		print(table,quote=FALSE,right=TRUE,na.print="",print.gap=2)

		cat("\n")
	}
}

