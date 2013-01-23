\name{mvmeta.control}
\alias{mvmeta.control}


\title{ Ancillary Parameters for Controlling the Fit in mvmeta Models }

\description{
This internal function sets the parameter options used for fitting meta-analytical models, commonly to pre-specified default values. It is usually internally called by \code{\link{mvmeta.fit}}.
}

\usage{
mvmeta.control(optim=list(), igls.niter=10)
}

\arguments{

  \item{optim }{ List of parameters passed to the \code{control} argument of the function \code{\link{optim}}, which perform the quasi-Newton optimization in likelihood-based random-effects models. See \code{\link{optim}}.}
  \item{igls.niter }{ Number of iteration of the iterative generalized least square algorithm to be run in the hybrid optimization procedure used in likelihood-based random-effects model.}
}

\details{
The control argument of \code{\link{mvmeta}} is by default passed to \code{\link{mvmeta.fit}}, which uses its elements as arguments of \code{mvmeta.control}.

The function will be eventually extended through other arguments related to different aspects of the fitting procedures.
}

\value{
A list with components named as the arguments.
}

\references{
Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/statmed2012.html}{\bold{here}}].
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\seealso{ 
See \code{\link{mvmeta}}. See also \code{\link{glm.control}}. See the different \code{\link[=mvmeta.ml]{fitting functions}} for details on each parameter. See \code{\link{mvmeta-package}} for an overview of this modelling framework.
}

\examples{
# RUN THE MODEL
model <- mvmeta(cbind(PD,AL)~pubyear,S=berkey98[5:7],data=berkey98)

# PRINT THE ITERATIONS (SEE ?optim)
model <- mvmeta(cbind(PD,AL)~pubyear,S=berkey98[5:7],data=berkey98,
  control=list(optim=list(trace=6,REPORT=1),igls.niter=20))
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
