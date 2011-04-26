\name{fitted.mvmeta}
\alias{fitted.mvmeta}


\title{ Model Fitted Values for mvmeta Objects }

\description{
This method function computes model fitted values for objects of class \code{"mvmeta"}.
}

\usage{
\method{fitted}{mvmeta}(object, format=c("matrix","list"), na.action, ...)
}

\arguments{

  \item{object}{ an object of class \code{"mvmeta"}.}
  \item{format}{ the format for the returned results. See Value.}
  \item{na.action}{ a function which indicates what should happen when the data contain \code{NA}s. The default to the value saved in the \code{"mvmeta"} object. See Note.}
  \item{\dots}{ further arguments passed to or from other methods.}

}

\value{
A \eqn{m} x \eqn{k} matrix or a \eqn{m}-length list of \eqn{k}-length vectors. In simple meta-analysis with no predictors, the fitted values are identical to model coefficients.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
This method function has a slightly different behaviour in the presence of missing values, if compared to the same method for other classes. Specifically, the fitted values are set to \code{NA} if at least one of the predictor values is missing in meta-regression models. In constrast, the function computes the fitted values in the presence of missing outcomes in \code{y} or entries of the within-study (co)variance matrices \code{S} in \code{\link{mvmeta}}.
}

\seealso{
The function \code{\link[mvmeta:predict.mvmeta]{predict}} computes the predicted values, optionally for different values of the predictors in meta-regression models. The function \code{\link[mvmeta:blup.mvmeta]{blup}} provides best linear unbiased pridictions for the set of studies included in the analysis. See \code{\link{mvmeta-package}} and \code{\link{mvmeta}} for an overview of the package and modelling framework.
}

\examples{
# RUN THE MODEL 
model <- mvmeta(y=berkey98[4:5],S=berkey98[6:8])

# FITTED VALUES
fitval <- fitted(model,format="list")
fitval
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
