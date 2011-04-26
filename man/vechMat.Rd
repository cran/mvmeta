\name{vechMat}
\alias{vechMat}
\alias{xpndMat}
\alias{expMat}
\alias{logMat}


\title{ Internal Functions for Package mvmeta }

\description{
These functions are used to apply some specific matrix operations. They are called internally by other functions within the \pkg{mvmeta} package, and are not meant to be directly run by the users (see Note below).
}

\usage{

vechMat(mat, diag=TRUE)

xpndMat(vech)

expMat(x)

logMat(x, positive=FALSE)

}

\arguments{

  \item{mat }{ a matrix.}
  \item{vech }{ a vector.}
  \item{x }{ a matrix or vector.}
  \item{diag }{ a logical switch indicating if the diagonal entries must be included.}
  \item{positive }{ a logical switch indicating if positive-definiteness must be forced.}

}

\details{
\code{vechMat} returns a vector with the lower triangle entries of the matrix, taken by column. \code{xpndMat} applied the reverse transformation.

\code{expMat} and \code{logMat} compute the matrix logarithm or exponential, applying these operations through eigenvalues decomposition of the matrix. For log-transformation, if \code{positive=TRUE}, the eigenvalues smaller than \code{1e-09} are set to this value.
}

\value{
A matrix or vector.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{As stated earlier, these functions are called internally within the \pkg{mvmeta} package, and are not meant to be used directly. In particular, their code does not contain any check on the arguments provided, which are expected in specific formats. Although not exported (included in the namespace), they are documented in order to facilitate inspection and intepretation of the package code. Their code is reported below.
}

\seealso{ 
See \code{\link{mvmeta-package}} and \code{\link{mvmeta}} for an overview of the package and modelling framework.
}

\examples{
vechMat <- function(mat, diag=TRUE) {
	if(!is.matrix(mat)) mat <- as.matrix(mat)
	return(mat[lower.tri(mat,diag=diag)])
}

xpndMat <- function(vech) {
	dim <- (-1 + sqrt(1 + 8 * length(vech)))/2
	mat <- matrix(nrow=dim,ncol=dim)
	mat[lower.tri(mat,diag=TRUE)] <- as.matrix(vech)
	mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
	return(mat)
}

expMat <- function(x) {
	if(is.data.frame(x)) x <- as.matrix(x)
	if(is.matrix(x)) {
		logA <- x
	} else logA <- xpndMat(x)
	eig <- eigen(logA)
	B <- diag(exp(eig$val),nrow(logA))
	A <- eig$vec%*%tcrossprod(B,eig$vec)
	return(A)
}

logMat <- function(x, positive=FALSE) {
	if(is.data.frame(x)) x <- as.matrix(x)
	if(is.matrix(x)) {
		expA <- x
	} else expA <- xpndMat(x)
	eig <- eigen(expA)
	if(positive) eig$val[eig$val<10^-9] <- 10^-9
	B <- diag(log(eig$val),nrow(expA))
	A <- eig$vec%*%tcrossprod(B,eig$vec)
	return(A)
}

}

\keyword{internal}
\keyword{manip}
