\name{gradchol}
\alias{gradchol}


\title{ Internal Functions for Package mvmeta }

\description{
This functions computes the partial derivatives of a symmetric matrix with respect to the entries of its Cholesky decomposition.
}

\usage{
gradchol(i, k, U)
}

\arguments{

  \item{i }{ entry position in the vector obtained by the upper triangle entries of the Cholesky decomposition of the matrix, taken by row.}
  \item{k }{ dimension of the symmetric matrix.}
  \item{U }{ upper triangular matrix obtained through Cholesky decomposition.}

}

\details{
#.
}

\value{
A symmetrix matrix representing the partial derivative of the original matrix with respect to the entry \code{i}.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
As stated earlier, these functions are called internally within the \pkg{mvmeta} package, and are not meant to be used directly. In particular, their code does not contain any check on the arguments provided, which are expected in specific formats. Although not exported (included in the namespace), they are documented in order to facilitate inspection and intepretation of the package code. Their code is reported below.
}

\seealso{ 
See \code{\link{mvmeta-package}} and \code{\link{mvmeta}} for an overview of the package and modelling framework.
}

\examples{
gradchol <- function(i,k,U) {
	ind1 <- rep(1:k,k:1)
	ind2 <- unlist(sapply(1:k,seq,to=k))
	C <- D <- E <- diag(0,k)
	C[ind2[i],] <- D[,ind2[i]] <- U[ind1[i],]
	E[ind2[i],] <- E[,ind2[i]] <- 1
	F <- E*C+E*D
	return(F)
}

}

\keyword{internal}
\keyword{manip}
