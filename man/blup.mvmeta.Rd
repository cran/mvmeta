\name{blup.mvmeta}
\alias{blup.mvmeta}
\alias{blup}

\title{ Best Linear Unbiased Prediction for mvmeta Objects }

\description{
This method function computes model best linear unbiased prediction for objects of class \code{"mvmeta"}. The function also returns standard errors, prediction intervals or the entire (co)variance matrix of the predicted outcome parameters. The original generic method is provided by the function \code{\link[metafor:blup]{blup}} within the package \pkg{metafor}.
}

\usage{
\method{blup}{mvmeta}(object, se=FALSE, pi=FALSE, vcov=FALSE, pi.level=0.95,
	format=c("matrix","list"), aggregate=c("stat","y"), na.action, ...)
}

\arguments{

  \item{object }{ an object of class \code{"mvmeta"}.}
  \item{se }{ logical switch indicating if standard errors must be included.}
  \item{pi }{ logical switch indicating if prediction intervals must be included.}
  \item{vcov }{ logical switch indicating if the (co)variance matrix must be included.}
  \item{pi.level }{ a numerical value between 0 and 1, specifying the confidence level for the computation of prediction intervals.}
  \item{format }{ the format for the returned results. See Value.}
  \item{aggregate }{ when \code{format="matrix"} and \code{se} or \code{ci} are required, the results may be aggregated by statistic or by outcome. See Value.}
  \item{na.action }{ a function which indicates what should happen when the data contain \code{NA}s. The default to the value saved in \code{object}. See Note.}
  \item{\dots }{ further arguments passed to or from other methods.}

}

\details{
For random-effects models, predictions are given by the sum of estimated outcome parameters from the fixed part of the model, plus study-specific deviations predicted as random effects given the between-study outcome parameters distribution.

Predicted outcome parameters from \code{blup} are a shrunk version of study-specific outcome realizations, where study-specific estimates borrow strength from the assumption of an underlying multivariate distribution of study parameters. In practice, the results from \code{blup} represent a weighted average between population outcome parameters (estimated by the fixed part of the model) and study-specific estimates. The weights depend from the relative size of the within and between-study covariance matrices \eqn{S} and \eqn{Psi}.

In the presence of missing values in the study-specific estimated outcome parameters \eqn{y} or in the related (co)variance matrix \eqn{S} of the fitted model, correspondent values of point estimates and covariance terms are set to 0, while the variance terms are set to \code{1e+10}. In this case, in practice, the study-specific estimates do not provide any information (their weight is virtually 0), and the prediction tends to the value returned by \code{\link[mvmeta:predict.mvmeta]{predict}} with \code{interval="prediction"}. See also Note below.

Fixed-effects models do not assume study-specific random effects, and the results of \code{blup} for these models are identical to \code{\link[mvmeta:predict.mvmeta]{predict}} with \code{interval="confidence"}.
}

\value{
The results may be aggregated in matrices (the default), or returned as lists, depending on the argument \code{format}. For multivariate models, the aggregation is ruled by the argument \code{aggregate}, and the results may be grouped by statistic or by outcome. If \code{vcov=TRUE}, lists are always returned.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
Similarly to \code{\link[mvmeta:predict.mvmeta]{predict}} or \code{\link[mvmeta:fitted.mvmeta]{fitted}}, the \code{blup} values are set to \code{NA} if at least one of the predictor values is missing in meta-regression models. Also, the function computes the predicted values in the presence of missing outcomes parameters estimates or entries of the within-study (co)variance matrices (left-hand side of \code{formula} and \code{S} in \code{\link{mvmeta}}, respectively). Interestingly, the results may be different from \code{\link[mvmeta:predict.mvmeta]{predict}} if at least some of the outcomes parameters are included in estimation, as the prediction also depends on the random part of the model. Specifically, the function uses information from the between-study covariance to predict missing outcome parametes given the observed ones.
}

\seealso{
See also the generic method function \code{\link[metafor:blup]{blup}}. The method function \code{\link[mvmeta:predict.mvmeta]{predict}} provides model predictions, optionally for a new set of predictor values in meta-regression models, while \code{\link[mvmeta:fitted.mvmeta]{fitted}} computes the fitted values. See \code{\link{mvmeta-package}} and \code{\link{mvmeta}} for an overview of the package and modelling framework.
}

\examples{
# RUN THE MODEL
model <- mvmeta(cbind(PD,AL)~pubyear,S=berkey98[6:8],data=berkey98)

# BLUP AND SE
blup(model,se=TRUE)
# SAME AS ABOVE, AGGREGATED BY OUTCOME
blup(model,se=TRUE,aggregate="y")
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}