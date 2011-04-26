\name{summary.mvmeta}
\alias{print.mvmeta}
\alias{summary.mvmeta}

\title{ Summary for mvmeta Objects }

\description{
Print and summary method functions for objects of class \code{"mvmeta"}.
}

\usage{
\method{print}{mvmeta}(x, ci.level=0.95, digits=4, ...)

\method{summary}{mvmeta}(object, ci.level=0.95, digits=4, ...)
}

\arguments{
  \item{x, object }{ an object of class \code{"mvmeta"}.}
  \item{ci.level }{ a numerical value between 0 and 1, specifying the confidence level for the computation of confidence intervals.}
  \item{digits }{ an integer specyfing the number of digits to which printed results must be rounded.}
  \item{\dots }{ further arguments passed to or from other methods.}
}

\details{
The \code{print} method only returns the main results, basically the estimates (with standard errors, z and p-values and confidence intervals for the fixed-effects \eqn{beta}), plus the components of the between-study (co)variance matrix \eqn{Psi} for random-effects models, namely standard deviations (square root of the diagonal) and correlations.

The \code{summary} method also prints information on the number of studies and outcome parameters included, the estimation method, number of observations, fixed and random parameters and fit statistics (log-likelihood, AIC, BIC).

Significance legend for the printed results: \code{'***'} < 0.001, \code{'**'} < 0.01, \code{'*'} < 0.05, \code{'.'} < 0.10.
}

\value{
These functions only print results from the model, and do not return any object.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\seealso{
See \code{\link{mvmeta-package}} and \code{\link{mvmeta}} for an overview of the package and modelling framework.
}


\examples{
# RUN THE MODEL 
model <- mvmeta(y=berkey98[4:5],S=berkey98[6:8])

# SIMPLE PRINT
model
# DEFINE DIGITS
print(model,digit=2)
# SUMMARY WITH 80TH CONFIDENCE INTERVALS
summary(model,ci.level=0.80)
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
