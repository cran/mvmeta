\name{kXlistmk}
\alias{kXlistmk}


\title{ Internal Functions for Package mvmeta }

\description{
This function creates the design matrix for the fixed part of the multivariate meta-analytic model.
}

\usage{
kXlistmk(X, cen, nalist, m, k)
}

\arguments{

  \item{X }{ a \eqn{m} x \eqn{p-1} matrix storing the predictor values for each of the \eqn{m} studies.}
  \item{cen }{ centering values for the predictors. Accepted formats are a \eqn{(p-1)}-length vector, or logical. If \code{TRUE}, centered on predictors means, if \code{FALSE} left uncentered.}
  \item{nalist }{ a list of \eqn{k}-length study-specific logical vectors, \code{TRUE} when the outcome parameter is not missing.}
  \item{m, k }{ numeric scalars: number of studies, number of outcome parameters.}

}

\details{
The function produces a list of desing matrix for multivariate meta-analysis and meta-regression models. In the absence of missing values, each matrix is produced by a Kronecker product between an identity matrix of dimension \eqn{k} and the transpose of the vector of study-level predictors. In the presence of missing outcome parameters, the related row of the identity matrix is excluded. For simple multivariate meta-analysis without predictors, the vector of study-level predictors is simply a numeric scalar of 1, representing the intercept. The final desing matrix is therefore an identity matrix itself.
}

\value{
A list of design matrices.
}

\references{
Berkey, C. S., D. C. Hoaglin, et al. (1998). "Meta-analysis of multiple outcomes by regression with random effects." \emph{Statistics in Medicine}. \bold{17}(22):2537-2550.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
As stated earlier, these functions are called internally within the \pkg{mvmeta} package, and are not meant to be used directly. In particular, their code does not contain any check on the arguments provided, which are expected in specific formats. Although not exported (included in the namespace), they are documented in order to facilitate inspection and intepretation of the package code. Their code is reported below.
}

\seealso{ 
See \code{\link{mvmeta-package}} and \code{\link{mvmeta}} for an overview of the package and modelling framework.
}

\examples{
kXlistmk <- function(X,cen,nalist,m,k) {
	# IF X DOES NOT EXIST, 
	if(is.null(X)) {
		kXlist <- mapply(function(i,na) diag(1,k)[na,,drop=FALSE],
			seq(m),nalist,SIMPLIFY=FALSE)
	# IF INSTEAD X IS A MATRIX
	} else {
		# CENTERING
		kXlist <- cbind(1,sweep(as.matrix(X),2,cen))
		# EXPAND
		kXlist <- mapply(function(i,na) 
			diag(1,k)[na,,drop=FALSE]\%x\%kXlist[i,,drop=FALSE],
			seq(m),nalist,SIMPLIFY=FALSE)
	}
	return(kXlist)
}

}

\keyword{internal}
\keyword{manip}
