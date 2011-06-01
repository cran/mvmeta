mvmeta <-
function(formula, S, data, method="reml", lab, contrasts, na.action, control) {

##########################################################################
# PREPARE THE DATA

	# CREATE THE CALL
	call  <- match.call()
	mcall <- match.call(expand.dots=FALSE)
	mn <- match(c("formula","data","constrasts"), names(mcall), 0L)
	mcall <- mcall[c(1L, mn)]
  	mcall$drop.unused.levels <- TRUE
	# HERE KEEP THE MISSING
  	mcall$na.action <- "na.pass"
	mcall[[1L]] <- as.name("model.frame")

	# CREATE FORMULA IF NOT PROVIDED (FOR SIMPLE META-ANALYSIS)
	if(missing(data)) {
		data <- sys.frame(sys.parent())
	} else if(!is.data.frame(data)) stop("'data' must be a data.frame")
	if(class(eval(mcall[[mn[1]]],data))!="formula") {
		call[[mn[1]]] <- mcall[[mn[1]]] <- formula <- 
			as.formula(paste(deparse(substitute(formula)),"~ 1"))
	}

	# CREATE y AND X FROM formula (FOLLOW lm FUNCTION)
	model <- eval(mcall,parent.frame())
	terms <- attr(model,"terms")
	# DESIGN MATRIX (NOT EXPANDED YET)
	if(missing(contrasts)) contrasts <- NULL
	X <- model.matrix(terms,model,contrasts)
	int <- any(colnames(X)=="(Intercept)")
	if(int) colnames(X)[colnames(X)=="(Intercept)"] <- "(Int)"
	# KEEP MATRIX FORMAT ALSO WITH UNIVARIATE RESPONSE
	y <- as.matrix(model.response(model,"numeric"))

	# NUMBER OF OUTCOMES AND PREDICTORS
	k <- ncol(y)
	p <- ncol(X)

	# CREATE S AND lab
	S <- eval(call$S,data)
	if(missing(lab)) {
		lab <- list()
	} else lab <- eval(call$lab,data)

	# FIRST CHECKS
	if(!is.null(model.offset(model))) stop("an offset is not allowed")
	if(is.empty.model(terms)) stop("an empty model is not allowed")
	if(is.null(model.response(model))) stop("response needed in formula")
	if(is.data.frame(S)) S <- as.matrix(S)
	if(missing(na.action)) na.action <- getOption("na.action")
	if(missing(control)) {
		control <- list()
	} else if(!is.list(control)) stop("'control' must be a list")
	mvmeta.check(y, S, X, method, lab, na.action)

	# LABELS
	mlab <- lab$mlab
	klab <- lab$klab
	plab <- lab$plab
	if(is.null(mlab)) {
		if(!is.null(rownames(y))) {
			mlab <- rownames(y)
		} else mlab <- paste("study",seq(nrow(y)),sep="")
	}
	rownames(y) <- rownames(X) <- mlab
	if(is.null(klab)) {
		labk <- paste("y",seq(k),sep="")
		if(is.null(colnames(y))) {
			colnames(y) <- labk
		} else colnames(y)[colnames(y)==""] <- labk[colnames(y)==""]
		klab <- colnames(y)
	} else colnames(y) <- klab
	if(is.null(plab)) {
		plab <- colnames(X)
	} else colnames(X) <- plab

	# TRANSFORM y IN LIST OF VECTORS
	y <- lapply(seq(nrow(as.matrix(y))), function(i) {
		yi <- as.matrix(y)[i,]
		return(yi)})

	# TRANFORM S IN A LIST OF MATRICES
	if(!is.list(S)) {
		if(length(dim(S))==3) {
			S <- lapply(seq(dim(S)[3]), function(i) {
				Si <- S[,,i]
				dimnames(Si) <- list(klab,klab)
				return(Si)})
		} else S <- lapply(seq(nrow(as.matrix(S))), function(i) {
			Si <- xpndMat(as.matrix(S)[i,])
			dimnames(Si) <- list(klab,klab)
			return(Si)})
	} else S <- lapply(S,function(x) {
		x <- as.matrix(x)
		dimnames(x) <- list(klab,klab)
		return(x)})
	names(S) <- mlab

	# MISSINGNESS: 
	#	CREATE A LOGICAL VECTOR OF STUDIES WITH NON-MISSING
	# 	CREATE A LIST WITH THE MISSING PATTERN (CONSIDERING y,S,X)
	nastudy <- sapply(y,function(x) !all(is.na(x)))
	nastudy[sapply(S,function(x) all(is.na(rowSums(x)+colSums(x))))] <- FALSE
	Xna <- !is.na(X)
	if(any(colSums(Xna)==0)) {
		stop("at least one predictor in 'X' is completely missing")
	}
	nastudy[rowSums(Xna)<ncol(X)] <- FALSE
	nalist <- mapply(function(y,S) {
		na <- !is.na(y)
		na[unique(which(is.na(S),arr.ind=TRUE))] <- FALSE
		return(na)},y[nastudy],S[nastudy],SIMPLIFY=FALSE)
	if(sum(nastudy)<2) stop("at least 2 non-missing studies are required")
	if(na.action=="na.fail" && !all(c(unlist(nalist),nastudy))) {
		stop("missing values in 'y' or 'X'")
	}

	# CREATE y AND S LISTS ELIMINATING MISSING
	ylist <- mapply(function(y,na) y[na],y[nastudy],nalist,SIMPLIFY=FALSE)
	Slist <- mapply(function(S,na) S[na,na,drop=FALSE],
		S[nastudy],nalist,SIMPLIFY=FALSE)

	# NUMBER OF STUDIES AND OBSERVATIONS
	m <- sum(nastudy)
	nobs <- sum(unlist(nalist))
	
	# CHECK X
	if(rankMatrix(X[nastudy,,drop=FALSE])<p) {
		stop("rank-deficient design matrix: check predictors and missing pattern",sep="")
	}
	# CREATE kXlist
	kXlist <- mapply(function(i,na) {diag(1,k)[na,,drop=FALSE]%x%
		X[nastudy,,drop=FALSE][i,,drop=FALSE]},
		seq(m),nalist,SIMPLIFY=FALSE)

##########################################################################
# FIXED-EFFECTS 

	if(method=="fixed") {

		# BETWEEN-STUDY (CO)VARIANCE MATRIX SET TO 0
		Psi <- diag(0,k)

		# OBTAIN beta BY GLS
		Sigmalist <- mapply(function(S,na) S+Psi[na,na,drop=FALSE],
			Slist,nalist,SIMPLIFY=FALSE)
		Ulist <- lapply(Sigmalist,chol)
		invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
		invtUXlist <- mapply(function(invU,kX) crossprod(invU,kX),
			invUlist,kXlist,SIMPLIFY=FALSE)
		invtUylist <- mapply(function(invU,y) crossprod(invU,y),
			invUlist,ylist,SIMPLIFY=FALSE)
		invtUX <- rbindlist(invtUXlist)
		invtUy <- rbindlist(invtUylist)
		beta <- as.numeric(qr.solve(invtUX,invtUy))

		# COMPUTE (CO)VARIANCE MATRIX OF beta
		qrinvtUX <- qr(invtUX)
		R <- qr.R(qrinvtUX)
		Qty <- qr.qty(qrinvtUX,invtUy)
		vcov <- tcrossprod(backsolve(R,diag(1,ncol(invtUX))))

		# OTHER OBJECTS
		logLik <- as.numeric(-0.5*(nobs*log(2*pi) +
			crossprod(invtUy-invtUX%*%beta))+
			-sum(sapply(Ulist,function(U) sum(log(diag(U))))))
		nranpar <- 0
		convergence <- NULL
	}

##########################################################################
# LIKELIHOOD METHODS (REML AND ML)

	if(method%in%c("reml","ml")) {

		# PRODUCE INITIAL VALUES
		Psi <- diag(0.001,k)
		niter <- 10
		for(i in 1:niter) {
			Psi <- mvmeta.igls(Psi,ylist,Slist,kXlist,nalist,k,m)
		}

		# GENERATE par AND RUN THE MODEL
		par <- vechMat(t(chol(Psi)))
		obj <- as.name(paste("mvmeta",method,sep="."))
		fgrad <- as.name(paste("mvmeta",method,"grad",sep="."))
		# FORCE TO MAXIMIZATION
		control$fnscale <- -1

		fit <- optim(par,eval(obj),eval(fgrad),ylist=ylist,Slist=Slist,
			kXlist=kXlist,nalist=nalist,nobs=nobs,k=k,method="BFGS",
			control=control)

		# Psi: ESTIMATED BETWEEN-STUDY (CO)VARIANCE MATRIX
		Psi <- matrix(0,k,k)
		Psi[lower.tri(Psi,diag=TRUE)] <- fit$par
		Psi <- tcrossprod(Psi)
	
		# OBTAIN beta BY GLS
		Sigmalist <- mapply(function(S,na) S+Psi[na,na,drop=FALSE],
			Slist,nalist,SIMPLIFY=FALSE)
		Ulist <- lapply(Sigmalist,chol)
		invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
		invtUXlist <- mapply(function(invU,kX) crossprod(invU,kX),
			invUlist,kXlist,SIMPLIFY=FALSE)
		invtUylist <- mapply(function(invU,y) crossprod(invU,y),
			invUlist,ylist,SIMPLIFY=FALSE)
		invtUX <- rbindlist(invtUXlist)
		invtUy <- rbindlist(invtUylist)
		beta <- as.numeric(qr.solve(invtUX,invtUy))

		# COMPUTE (CO)VARIANCE MATRIX OF beta
		qrinvtUX <- qr(invtUX)
		R <- qr.R(qrinvtUX)
		Qty <- qr.qty(qrinvtUX,invtUy)
		vcov <- tcrossprod(backsolve(R,diag(1,ncol(invtUX))))

		# LOGLIKELIHOOD AND CONVERGENCE
		logLik <- fit$value
		nranpar <- length(fit$par)
		convergence <- fit$convergence
}

##########################################################################

	# DF AND LABELS
	df <- list(nobs=nobs,fixed=length(beta),random=nranpar)
	lab <- list(mlab=mlab,klab=klab,plab=plab)
	dim <- list(m=m,k=k,p=p)
	
	# PUT NAMES ON beta, vcov ANS Psi
	plablong <- paste(rep(klab,each=length(plab)),rep(plab,k),sep=".")
	names(beta) <- plablong
	dimnames(vcov) <- list(plablong,plablong)
	dimnames(Psi) <- list(klab,klab)	

	results <- list(beta=beta,vcov=vcov,Psi=Psi,method=method,
		logLik=logLik,df=df,dim=dim,lab=lab,int=int,call=call,model=model,
		S=S,na.action=na.action,contrasts=contrasts,
		convergence=convergence)

	class(results) <- "mvmeta"
	return(results)
}

