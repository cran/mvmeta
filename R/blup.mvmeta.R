blup.mvmeta <-
function(object, se=FALSE, pi=FALSE, vcov=FALSE, pi.level=0.95,
	format=c("matrix","list"), aggregate=c("stat","y"), na.action, ...) {

	# CHECKS
	if(pi.level<=0||pi.level>=1) stop("'pi.level' must be within 0 and 1")
	format <- match.arg(format,c("matrix","list"))
	aggregate <- match.arg(aggregate,c("stat","y"))
	if(missing(na.action)) na.action <- object$na.action
	na.action <- match.arg(na.action,c("na.omit","na.exclude",
		"na.fail","na.pass"))

#########################################################################

	# RE-CREATE OBJECTS
	# HERE nalist2 INCLUDES ALL THE OBS, ALSO MISSING
	# YLIST INCLUDING 0 WHEN NA
	# SLIST INCLUDING 0 AND 10^12 IN COVARIANCES AND VARIANCES,
	#	RESPECTIVELY, MATCHING THE PATTERN IN y
	# kXlist INCLUDES ALL THE OBS, ALSO MISSING
	model <- object$model
	y <- as.matrix(model.response(model,"numeric"))
	y <- lapply(seq(nrow(as.matrix(y))), function(i) as.matrix(y)[i,])
	terms <- terms(model)
	X <- model.matrix(terms,model,object$contrast)
	S <- object$S
	nalist <- mapply(function(y,S) {
		na <- !is.na(y)
		na[unique(which(is.na(S),arr.ind=TRUE))] <- FALSE
		return(na)},y,S,SIMPLIFY=FALSE)
	nalist2 <- lapply(nalist,function(x) rep(TRUE,length(x)))
	ylist <- mapply(function(y,na) {
		y[!na] <- 0
		return(y)},y,nalist,SIMPLIFY=FALSE)
	Slist <- mapply(function(S,na) {
		S1 <- diag(10^10,length(na))
		S1[na,na] <- S[na,na]
		return(S1)},S,nalist,SIMPLIFY=FALSE)
	kXlist <- mapply(function(i,na) {diag(1,object$dim$k)[na,,drop=FALSE]%x%
		X[i,,drop=FALSE]},seq(length(ylist)),nalist2,SIMPLIFY=FALSE)

		
	# PREPARE ADDITIONAL OBJECTS
	Sigmalist <- lapply(Slist, function(S) S+object$Psi)
	Ulist <- lapply(Sigmalist,chol)
	invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
	reslist <- mapply(function(y,kX) res <- y-kX%*%object$beta,
		ylist,kXlist,SIMPLIFY=FALSE)
	Psi <- object$Psi

	# COMPUTE BLUP
	bluplist <- mapply(function(invU,res,kX) {
		blup <- as.numeric(kX%*%object$beta + 
			Psi%*%(tcrossprod(invU))%*%res)
		names(blup) <- object$lab$klab
		return(blup)},invUlist,reslist,kXlist,SIMPLIFY=FALSE)

	# COMPUTE VCOV (SUM OF COMPONENTS), SE AND CI
	zvalpi <- qnorm((1-pi.level)/2,lower.tail=FALSE)
	vcovlist <- mapply(function(y,S,kX,invU) {
		vcov <- kX%*%tcrossprod(object$vcov,kX) +
			Psi - Psi%*%(tcrossprod(invU)%*%Psi)
		dimnames(vcov) <- list(object$lab$klab,object$lab$klab)
		return(vcov)},ylist,Slist,kXlist,invUlist,SIMPLIFY=FALSE)
	selist <- lapply(vcovlist, function(x) sqrt(diag(x)))
	pilblist <- mapply(function(blup,se) blup-zvalpi*se,
			bluplist,selist,SIMPLIFY=FALSE)
	piublist <- mapply(function(blup,se) blup+zvalpi*se,
			bluplist,selist,SIMPLIFY=FALSE)

	# HANDLE MISSING
	Xna <- rowSums(!is.na(X))==object$dim$p
	if(na.action=="na.fail" && !all(Xna)) {
		stop("missing values in 'X'")
	}
	if(na.action=="na.omit") {
		bluplist <- bluplist[Xna]
		selist <- selist[Xna]
		pilblist <- pilblist[Xna]
		piublist <- piublist[Xna]
		vcovlist <- vcovlist[Xna]
		mlab <- object$lab$mlab[Xna]
	} else mlab <- object$lab$mlab

#########################################################################

	# IF VCOV WHEN MORE THAN 1 OUTCOME, SWITCH TO LIST
	if(vcov && object$dim$k>1) format <- "list"
	# IF ONLY 1 OUTCOME, FORCE THE AGGREGATION ON IT
	if(format=="matrix" && object$dim$k==1) {
		blup <- rbindlist(bluplist)
		if(se) blup <- cbind(blup,rbindlist(selist))
		if(pi) blup <- cbind(blup,rbindlist(pilblist),rbindlist(piublist))
		if(vcov) blup <- cbind(blup,rbindlist(vcovlist))
		dimnames(blup) <- list(mlab,c("blup","se","pi.lb",
			"pi.ub","vcov")[c(TRUE,se,pi,pi,vcov)])
	# IF NO OTHER STAT, FORCE THE AGGREGATION ON PREDICTION
	}else if(format=="matrix" && !se && !pi) {
		blup <- rbindlist(bluplist)
		rownames(blup) <- mlab
	# WHEN AGGREGATE ON STAT
	} else if (format=="matrix" && aggregate=="stat") {
		blup <- list(blup=rbindlist(bluplist))
		rownames(blup$blup) <- mlab
		if(se) {
			blup$se <- rbindlist(selist)
			rownames(blup$se) <- mlab
		}
		if(pi) {
			blup$pi.lb <- rbindlist(pilblist)
			rownames(blup$pi.lb) <- mlab
			blup$pi.ub <- rbindlist(piublist)
			rownames(blup$pi.ub) <- mlab
		}
	# WHEN AGGREGATE ON OUTCOME
	} else if (format=="matrix" && aggregate=="y") {
		blup <- list()
		for(j in seq(object$dim$k)) {
			templist <- lapply(seq(bluplist),function(i) {
				blup <- bluplist[[i]][j]
				if(se) blup <- cbind(blup,selist[[i]][j])
				if(pi) blup <- cbind(blup,pilblist[[i]][j],
					piublist[[i]][j])
				return(as.matrix(blup))})
			temp <- rbindlist(templist)
			dimnames(temp) <- list(mlab,
				c("blup","se","pi.lb","pi.ub")[c(TRUE,se,pi,pi)])
			blup[[j]] <- temp
		}
		names(blup) <- object$lab$klab
	# ALL THE OTHER COMBINATIONS, PRODUCE A LIST
	} else {
		blup <- mapply(function(bluparg,searg,pilbarg,piubarg,vcovarg) {
			temp <- list(blup=bluparg)
			if(se) temp$se <- searg
			if(pi) temp$pi.lb <- pilbarg
			if(pi) temp$pi.ub <- piubarg
			if(vcov) temp$vcov <- vcovarg
			return(temp)},
			bluplist,selist,pilblist,piublist,vcovlist,SIMPLIFY=FALSE)
		names(blup) <- mlab
	}

	return(blup)
}

