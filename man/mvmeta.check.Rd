\name{mvmeta.check}
\alias{mvmeta.check}


\title{ Internal Functions for Package mvmeta }

\description{
This function performs some checks on the arguments of the main function \code{\link{mvmeta}}.
}

\usage{
mvmeta.check(y, S, X, method, lab, cen, na.action)
}

\arguments{

See also the related arguments in \code{\link{mvmeta}}:

  \item{y }{ series of \eqn{k} estimated outcome parameters for each one of the \eqn{m} studies.}
  \item{S }{ series of within-study (co)variance matrices of estimated outcome parameters for each one of the \eqn{m} studies.}
  \item{X }{ series of \eqn{p-1} study-level predictors (excluding intercept) for each one of the \eqn{m} studies.}
  \item{method }{ estimation method.}
  \item{lab }{ list of label vectors.}
  \item{cen }{ centering values for the predictors in meta-regression.}
  \item{na.action }{ a function which indicates what should happen when the data contain \code{NA}s.}

}

\details{
The function run checks on format and dimensionality of each argument.
}

\value{
No result is produced, but an error message is returned if some argument does not match the expected format or dimension.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}
\note{
As stated earlier, these functions are called internally within the \pkg{mvmeta} package, and are not meant to be used directly. In particular, their code does not contain any check on the arguments provided, which are expected in specific formats. Although not exported (included in the namespace), they are documented in order to facilitate inspection and intepretation of the package code. Their code is reported below.
}

\seealso{ 
See \code{\link{mvmeta-package}} and \code{\link{mvmeta}} for an overview of the package and modelling framework.
}

\examples{
mvmeta.check <- function(y, S, X, method, lab, cen, na.action) {
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
	if(!method\%in\%c("fixed","ml","reml")) {
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

}

\keyword{internal}
