mvmeta.check <-
function(y, S, X, method, lab, cen, na.action) {

	# y, S AND X TYPES
	if(!is.numeric(y) && !is.list(y)) {
		stop("the argument 'y' must be a matrix, vector, list or dataframe")
	}
	if(!is.numeric(S) && !is.list(S) && !is.array(S)) {
		stop("the argument 'S' must be a matrix, vector, list, dataframe or array")
	}
	if(!is.null(X) && !is.numeric(X) && !is.list(X)) {
		stop("the argument 'y' must be a matrix, vector, list or dataframe")
	}

	# DIMENSIONS
	if(is.numeric(y)) {
		m <- nrow(as.matrix(y))
		k <- ncol(as.matrix(y))
	} else {
		m <- length(y)
		k <- nrow(as.matrix(y[[1]]))
	}
	if(m<2) stop("at least 2 studies are required")

	# CHECK ON y
	if(!is.numeric(y) && any(sapply(y,length)!=k)) {
			stop("studies have different number of outcomes")
	}
	# CHECK ON S
	if(is.numeric(S) && !is.matrix(S)) {
		if(k>1) stop("the argument S must provide a kxk covar matrix for m studies")
		if(length(S)!=m) stop("Dimension of 'y' must be the same of 'S'")
	} else if(is.matrix(S)) {
		if(ncol(S)!=k*(k+1)/2) stop("incorrect dimensions of argument 'S'")
		if(nrow(S)!=m) stop("Dimension of 'y' must be the same of 'S'")
	} else if(is.list(S)) {
		if(any(sapply(S,function(x) dim(as.matrix(x))!=k))) {
			stop("incorrect dimensions of argument 'S'")
		}
		if(length(S)!=m) stop("Dimension of 'y' must be the same of 'S'")
	} else {
		if(dim(S)!=c(k,k,m)) stop("incorrect dimensions of argument 'S'")
	}
	# CHECK ON X
	if(!is.null(X) && is.numeric(X)) {
		p <- 1+ncol(as.matrix(X))
		if(nrow(as.matrix(X))!=m) stop("incorrect dimensions of argument 'X'")
	} else if(!is.null(X)) {
		if(length(X)!=m) stop("incorrect dimensions of argument 'X'")
		p <- length(X[[1]])
		if(any(sapply(X,length)!=p)) stop("incorrect dimensions of argument 'X'")
	} else p <- 1

	# METHOD
	if(!method%in%c("fixed","ml","reml")) {
		stop("Methods allowed: 'fixed', 'ml', 'reml'")
	}

	# LAB
	if(!is.list(lab)) stop("argument 'lab' must be a list")
	if(!is.null(lab$mlab) && length(lab$mlab)!=m) {
		 stop("incorrect dimensions of argument 'mlab'")
	}
	if(!is.null(lab$klab) && length(lab$klab)!=k) {
		 stop("incorrect dimensions of argument 'klab'")
	}
	if(!is.null(lab$plab) && length(lab$plab)!=p-1) {
		 stop("incorrect dimensions of argument 'plab'")
	}

	#CEN
	if(!is.logical(cen) && !is.numeric(cen)) {
		stop("argument 'cen' must be logical or a numeric vector")
	}
	if(!is.null(X) && is.numeric(cen)) {
		if(length(cen)!=p-1) stop("incorrect dimensions of argument 'cen'")
	}

	# NA.ACTION
		format <- match.arg(na.action,c("na.omit","na.exclude",
		"na.fail","na.pass"))
}

