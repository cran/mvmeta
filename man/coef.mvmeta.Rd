\name{coef.mvmeta}
\alias{coef.mvmeta}
\alias{vcov.mvmeta}


\title{ Model Coefficients and their (Co)Variance Matrix for mvmeta Objects }

\description{
These method functions extract the estimated model coefficients and their (co)variance matrix for objects of class \code{"mvmeta"}.
}

\usage{

\method{coef}{mvmeta}(object, ...)

\method{vcov}{mvmeta}(object, ...)

}

\arguments{
  \item{object}{ an object of class \code{"mvmeta"}.}
  \item{\dots}{ further arguments passed to or from other methods.}
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\seealso{
See \code{\link{mvmeta-package}} and \code{\link{mvmeta}} for an overview of the package and modelling framework.
}

\examples{
# RUN THE MODEL 
model <- mvmeta(cbind(PD,AL)~pubyear,S=berkey98[6:8],data=berkey98)

# EXTRACT COEFFICIENTS
coef(model)
# EXTRACT (CO)VARIANCE MATRIX
vcov(model)
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
