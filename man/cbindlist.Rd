\name{cbindlist}
\alias{cbindlist}
\alias{rbindlist}
\alias{sumlist}
\alias{meanlist}


\title{ Internal Functions for Package mvmeta }

\description{
These functions are used to transform lists of matrices in a single matrix, applying specific operations. They are called internally by other functions within the \pkg{mvmeta} package, and are not meant to be directly run by the users (see Note below).
}

\usage{

cbindlist(list)

rbindlist(list)

sumlist(list)

meanlist(list)

}

\arguments{

  \item{list }{ A list of matrices with specific dimensions.}

}

\details{
\code{cbindlist} and \code{rbindlist} simply append by columns or rows the set of matrices. All the matrices are expected to have the same number of columns or rows, respectively.

\code{sumlist} and \code{meanlist} take the  entrywise sum or mean. All the matrices are expected to have the same dimensions.
}

\value{
A matrix.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
As stated earlier, these functions are called internally within the \pkg{mvmeta} package, and are not meant to be used directly. In particular, their code does not contain any check on the arguments provided, which are expected in specific formats. Although not exported (included in the namespace), they are documented in order to facilitate inspection and intepretation of the package code. Their code is reported below.
}

\seealso{ 
See \code{\link{mvmeta-package}} and \code{\link{mvmeta}} for an overview of the package and modelling framework.
}

\examples{
cbindlist <- function(list) {
	n <- length(list)
	res <- NULL
	for (i in seq(n)) res <- cbind(res, list[[i]])
	return(res)
}

rbindlist <-
function(list) {
	n <- length(list)
	res <- NULL
	for (i in seq(n)) res <- rbind(res, list[[i]])
	return(res)
}

sumlist <- function(list) {
	n <- length(list)
	res <- 0
	for (i in seq(n)) res <- res+list[[i]]
	return(res)
}

meanlist <- function(list) {
	n <- length(list)
	res <- 0
	for (i in seq(n)) res <- res+list[[i]]
	return(res/n)
}

}

\keyword{internal}
\keyword{manip}
