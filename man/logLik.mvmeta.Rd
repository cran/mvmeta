\name{logLik.mvmeta}
\alias{logLik.mvmeta}
\alias{AIC.mvmeta}
\alias{BIC.mvmeta}


\title{ Model Log-Likelihood and Information Criteria for mvmeta Objects }

\description{
These method functions extract the (restricted) log-likelihood, Akaike (AIC) and Bayesian (BIC) information criteria for objects of class \code{"mvmeta"}.
}

\usage{

\method{logLik}{mvmeta}(object, ...)

\method{AIC}{mvmeta}(object, ..., k = 2)

\method{BIC}{mvmeta}(object, ...)

}

\arguments{
  \item{object}{ an object of class \code{"mvmeta"}.}
  \item{k}{ numeric scalar: the penalty per parameter. The default \code{k=2} is the classical AIC.}
  \item{\dots}{ further arguments passed to or from other methods.}
}

\value{
For \code{logLik}, an object of class \code{"logLik"} with attributes, providing the (restricted) log likelihood of the model. See \code{\link{logLik}} for additional information.

For \code{AIC} and \code{BIC}, a numeric value corresponding to the selected information criteria. See \code{\link{AIC}}, \code{\link{BIC}} and \code{\link{extractAIC}} for additional information.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
The attributes of the object returned by \code{logLik} may display different information than that stored in the \code{"mvmeta"} model object. In particular, \code{nall} defines the total number of observations used for estimation (equal to \code{df$nobs} in the model object), while here \code{nobs} returns the same number, but substracted by the number of estimated coefficients for models estimated by REML. The number of total estimated parameters \code{df} is given by the sum of number of coefficients (\code{fixed}) and (co)variance terms (\code{random}) in the model. This choice is made in order to be coherent with related methods functions for other classes in the computation of AIC and BIC.
}

\section{Warnings}{
Particular care is needed when using such statistics for model comparison. In particular, models estimated by REML do not specify a true likelihood, and statistical theory underlying the specification of AIC and BIC does not fully apply. This is especially relevant for REML models estimated with a different set of predictors in the fixed part.
}

\seealso{ 
See \code{\link{mvmeta-package}} and \code{\link{mvmeta}} for an overview of the package and modelling framework.
}

\examples{
# RUN THE MODEL 
model <- mvmeta(cbind(PD,AL)~pubyear,S=berkey98[6:8],data=berkey98)

# AIC and BIC
AIC(model)
BIC(model)
# LOG-LIKELIHOOD
ll <- logLik(model)
ll
attributes(ll)
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
