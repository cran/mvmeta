\name{mvmeta.control}
\alias{mvmeta.control}


\title{ Ancillary Parameters for Controlling the Fit in mvmeta Models }

\description{
This internal function sets the parameter options used for fitting meta-analytical models, commonly to pre-specified default values. It is usually called by \code{\link{mvmeta.fit}}.
}

\usage{
mvmeta.control(fnscale=-1, \dots)
}

\arguments{

  \item{fnscale }{ An overall scaling to likelihood functions applied during optimization. If negative, turns the problem into a maximization problem. This should not be changed.}
  \item{\dots }{ Additional parametes. See \code{\link{optim}}.}

}

\details{
The control argument of \code{\link{mvmeta}} is by default passed to the control argument of \code{\link{mvmeta.fit}}, which uses its elements as arguments to \code{glm.control}.

In the present version, the functions simply specifies the parameter values to be passed to \code{\link{optim}}, which is called internally by \code{\link{mvmeta.ml}} and \code{\link{mvmeta.reml}} for controlling the optimization procedure in random-effects models.

The function will be eventually extended through other arguments related to different aspects of the fitting procedures.
}

\value{
A list with components named as the arguments.
}

\references{
Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. Epub ahead of print (doi 10.1002/sim.5471). [Freely available \href{http://onlinelibrary.wiley.com/doi/10.1002/sim.5471/pdf}{\bold{here}}].
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\seealso{ 
See \code{\link{mvmeta}}. See also \code{\link{glm.control}}. See \code{\link{mvmeta-package}} for an overview of this modelling framework.
}

\examples{
# RUN THE MODEL
model <- mvmeta(cbind(PD,AL)~pubyear,S=berkey98[5:7],data=berkey98)

# PRINT THE ITERATIONS (SEE ?optim)
model <- mvmeta(cbind(PD,AL)~pubyear,S=berkey98[5:7],data=berkey98,
  control=list(trace=6,REPORT=1))
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
