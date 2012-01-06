\name{mvmeta.check}
\alias{mvmeta.check}

\title{ Internal Functions for Package mvmeta }

\description{
This function performs some checks on the arguments of the main function \code{\link{mvmeta}}.
}

\usage{
mvmeta.check(y, S, X, method, lab, na.action)
}

\arguments{

See also the related arguments in \code{\link{mvmeta}}:

  \item{y }{ \eqn{m} x \eqn{k} matrix of estimated outcome parameters for each one of the \eqn{m} studies.}
  \item{S }{ series of within-study (co)variance matrices of estimated outcome parameters for each one of the \eqn{m} studies, allowed in different formats.}
  \item{X }{ \eqn{m} x \eqn{p} matrix of study-level predictors (including intercept) for each one of the \eqn{m} studies.}
  \item{method }{ estimation method.}
  \item{lab }{ list of label vectors.}
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
mvmeta.check <- function(y, S, X, method, lab, na.action) {

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
	if(!is.null(lab$plab) && length(lab$plab)!=p) {
		 stop("incorrect dimensions of argument 'plab'")
	}

	# NA.ACTION
		format <- match.arg(na.action,c("na.omit","na.exclude",
		"na.fail","na.pass"))
}

}

\keyword{internal}
