mvmeta.check <-
function(y, S, X, method, lab, na.action) {

	# S TYPES
	if(!is.numeric(S) && !is.list(S) && !is.array(S)) {
		stop("the argument 'S' must be a matrix, vector, list, dataframe or array")
	}

	# DIMENSIONS
	m <- nrow(y)
	k <- ncol(y)
	p <- ncol(X)
	if(m<2) stop("at least 2 studies are required")

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
	if(!is.null(lab$plab) && length(lab$plab)!=p) {
		 stop("incorrect dimensions of argument 'plab'")
	}

	# NA.ACTION
		format <- match.arg(na.action,c("na.omit","na.exclude",
		"na.fail","na.pass"))
}

