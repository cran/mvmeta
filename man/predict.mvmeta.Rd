\name{predict.mvmeta}
\alias{predict.mvmeta}


\title{ Predicted Values for mvmeta Objects }

\description{
This method function computes model prediction for objects of class \code{"mvmeta"}, optionally for a new set of predictor values in meta-regression models. The function also returns standard errors, confidence intervals or the entire (co)variance matrix of the predicted outcome parameters.
}

\usage{
\method{predict}{mvmeta}(object, newdata, se=FALSE, ci=FALSE, vcov=FALSE,
	interval=c("confidence","prediction"), ci.level=0.95,
	format=c("matrix","list"), aggregate=c("stat","y"), na.action, ...)
}

\arguments{

  \item{object }{ an object of class \code{"mvmeta"}.}
  \item{newdata }{ An optional data frame in which to look for variables values with which to predict from meta-regression models.}
  \item{se }{ logical switch indicating if standard errors must be included.}
  \item{ci }{ logical switch indicating if confidence intervals must be included.}
  \item{vcov }{ logical switch indicating if the (co)variance matrix must be included.}
  \item{interval }{ type of prediction. See Details.}
  \item{ci.level }{ a numerical value between 0 and 1, specifying the confidence level for the computation of confidence intervals.}
  \item{format }{ the format for the returned results. See Value.}
  \item{aggregate }{ when \code{format="matrix"} and \code{se} or \code{ci} are required, the results may be aggregated by statistic or by outcome. See Value}
  \item{na.action }{ a function which indicates what should happen when the data contain \code{NA}s. The default to the value saved in \code{object}. See Note.}
  \item{\dots }{ further arguments passed to or from other methods.}

}

\details{
For both fixed and random-effects models, estimated predictions are only based on the fixed part of the model, ignoring study-specific deviations, differently from \code{\link[mvmeta:blup.mvmeta]{blup}}.

For multivariate meta-analysis, the results are identical to the estimated coefficients and those returned by \code{\link[mvmeta:fitted.mvmeta]{fitted}}, but a single value is returned. Interpretation of fixed and random-effects models simply follows. In meta-regression models, results of random-effects models are interpreted as the estimated average outcome parameters of a set of studies characterized by the same specified values of predictors, sampled from the same (hypothetical) population the original studies included in the analysis have been selected from. Interpretation of predictions from fixed-effects meta-regression models is not straightforward (or even meaningful).

For random-effects models, if \code{interval="confidence"} (the default), standard errors, confidence intervals and (co)variance matrix of the predicted values are computed only using the estimated (co)variance matrix of the coefficients \eqn{beta}, coherently with the interpretation given above. If \code{interval="prediction"}, the estimated between-study (co)variance matrix \eqn{Psi} is also added. In this case, interpretation of the uncertainty reflects that of a new single study sampled from the same population.

}

\value{
The results may be aggregated in matrices (the default), or returned as lists, depending on the argument \code{format}. For aggregated results of a single prediction or meta-analysis with no predictor, vectors are returned. For multivariate models, the aggregation is ruled by the argument \code{aggregate}, and the results may be grouped by statistic or by outcome. If \code{vcov=TRUE}, lists are always returned. 
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
This method function has a slightly different behaviour in the presence of missing values, if compared to the same method for other classes. Specifically, the predicted values are set to \code{NA} if at least one of the predictor values is missing in meta-regression models. In constrast, the function computes the predicted values in the presence of missing outcomes parameters estimates or entries of the within-study (co)variance matrices (left-hand side of \code{formula} and \code{S} in \code{\link{mvmeta}}, respectively).
}

\seealso{ 
The method function \code{\link[mvmeta:blup.mvmeta]{blup}} provides best linear unbiased pridictions for the set of studies included in the analysis. The method function \code{\link[mvmeta:fitted.mvmeta]{fitted}} computes the fitted values. See \code{\link{mvmeta-package}} and \code{\link{mvmeta}} for an overview of the package and modelling framework.
}

\examples{
# RUN THE MODEL
model <- mvmeta(cbind(PD,AL)~pubyear,S=berkey98[6:8],data=berkey98)

# PREDICTED AVERAGED OUTOCOMES AND SE FROM YEAR 1985 TO 1987
newdata <- data.frame(pubyear=1985:1987)
predict(model,newdata,se=TRUE)
# SAME AS ABOVE, AGGREGATED BY OUTCOME
predict(model,newdata,se=TRUE,aggregate="y")
# SAME AS ABOVE, WITH PREDICTION INTERVALS
predict(model,newdata,se=TRUE,aggregate="y",interval="prediction")
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
