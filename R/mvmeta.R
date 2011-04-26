mvmeta <-
function(y, S, X=NULL, method="reml", lab, cen=FALSE, na.action, ...) {

##########################################################################
# PREPARE THE DATA

	# FIRST CHECKS
	if(is.data.frame(y)) y <- as.matrix(y)
	if(is.data.frame(S)) S <- as.matrix(S)
	if(is.data.frame(X)) X <- as.matrix(X)
	if(missing(lab)) lab <- list()
	if(missing(na.action)) na.action <- getOption("na.action")
	mvmeta.check(y, S, X, method, lab, cen, na.action)

	# TRANSFORM y IN LIST
	if(!is.list(y)) {
		mlab <- rownames(y)
		klab <- colnames(y)
		y <- lapply(seq(nrow(as.matrix(y))), function(i) {
			yi <- as.matrix(y)[i,]
			names(yi) <- klab
			return(yi)})
		names(y) <- mlab
	}

	# TRANFORM X IN MATRIX, CENTER AND EXPAND + LABELS
	if(is.list(X)) {X <- rbindlist(X)
	} else if(!is.null(X)) X <- as.matrix(X)

	# LABELS: PUT NAMES OF FIRST OBS OR CREATE THEM
	mlab <- lab$mlab
	klab <- lab$klab
	plab <- lab$plab
	if(is.null(mlab)) {
		if(!is.null(names(y))) { mlab <- names(y)
		} else {
			mlab <- paste("study",seq(length(y)),sep="")
			names(y) <- mlab
		}
	}
	if(is.null(klab)) {
		if(!is.null(names(y[[1]]))) { klab <- names(y[[1]])
		} else {
			klab <- paste("y",seq(length(y[[1]])),sep="")
			for(i in seq(length(y))) names(y[[i]]) <- klab
		}
	}
	if(is.null(plab)&&!is.null(X)) {
		if(!is.null(colnames(X))) { plab <- colnames(X)
		} else plab <- paste("beta",seq(ncol(as.matrix(X))),sep="")
	}
	plab <- c("int",plab)
	if(!is.null(X)) {
		rownames(X) <- mlab
		colnames(X) <- plab[-1]
	}

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
	if(!is.null(X)) {
		Xna <- !is.na(X)
		if(any(colSums(Xna)==0)) {
			stop("at least one predictor in 'X' is completely missing")
		}
		nastudy[rowSums(Xna)<ncol(X)] <- FALSE
	}
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

	# NUMBER OF STUDIES, OUTCOMES AND PREDICTORS
	m <- sum(nastudy)
	k <- length(y[[1]])
	p <- length(plab)
	nobs <- sum(unlist(nalist))
	
	# CREATE kXlist
	if(!is.null(X)) {
		if(!is.numeric(cen)) {
			if(cen) { cen <- colMeans(as.matrix(X),na.rm=T)
			} else cen <- rep(0,ncol(as.matrix(X)))
		}
		names(cen) <- plab[-1]
	} else cen <- FALSE
	kXlist <- kXlistmk(X[nastudy,,drop=FALSE],cen,nalist,m,k)

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
			Psi <- mvmeta.igls(Psi,ylist,Slist,kXlist,nalist,k,m,p)
		}

		# GENERATE par AND RUN THE MODEL
		par <- vechMat(t(chol(Psi)))
		obj <- as.name(paste("mvmeta",method,sep="."))
		fgrad <- as.name(paste("mvmeta",method,"grad",sep="."))

		model <- optim(par,eval(obj),eval(fgrad),ylist=ylist,Slist=Slist,
			kXlist=kXlist,nalist=nalist,nobs=nobs,k=k,method="BFGS",...)

		# Psi: ESTIMATED BETWEEN-STUDY (CO)VARIANCE MATRIX
		Psi <- matrix(0,k,k)
		Psi[lower.tri(Psi,diag=TRUE)] <- model$par
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
		logLik <- -model$value
		nranpar <- length(model$par)
		convergence <- model$convergence
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

	results <- list(y=y,S=S,X=X,beta=beta,vcov=vcov,
		Psi=Psi,method=method,logLik=logLik,df=df,dim=dim,lab=lab,cen=cen,
		na.action=na.action,convergence=convergence)

	class(results) <- "mvmeta"
	return(results)
}

