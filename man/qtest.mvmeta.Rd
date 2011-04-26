\name{qtest.mvmeta}
\alias{qtest.mvmeta}


\title{ Multivariate Cochran Q Test of Heterogeneity for mvmeta objects }

\description{
This method function perform a Cochran Q test of (residual) heterogeneity on multivariate meta-analytic models saved in objects of class \code{"mvmeta"}. The generic method is provided by \code{\link{qtest}}.
}

\usage{
\method{qtest}{mvmeta}(object, ...)
}

\arguments{

  \item{object}{ an object of class \code{"mvmeta"}.}
  \item{\dots}{ further arguments passed to or from other methods.}

}

\details{
In this context, the test assesses the null hypothesis that the variability in the multivariate distribution of the estimated outcome parameters is explained only in terms of sampling variability, measured by the within-study (co)variance matrices \eqn{S}. This is equal to test the hypothesis that the between-study covariance matrix \eqn{Psi} is a zero matrix, and there is no random deviation in study-specific estimates. In practice, the function re-runs a fixed-effects model and tests if the distribution of residuals is coherent with the null hypothesis.
}

\value{
A list object of class \code{"qtest.mvmeta"} with the following components:
  \item{Q }{ the test statistic, distributed under the null hypothesis as a Chi-square with degrees of freedom \code{df}.}
  \item{df }{ the degrees of freedom of the null distribution, equal to the number of observations used for estimation minus the number of coefficients in the fixed part of the model.}
  \item{pvalue }{ the p-value of the test.}
  \item{residual }{ logical switch indicating if a meta-regression model is assessed, meaning that the tested heterogeneity is residual.}
  \item{k }{ dimensionality of the test, that is the number of outcome parameters in the model.}
}

\references{
Berkey, C. S., D. C. Hoaglin, et al. (1998). "Meta-analysis of multiple outcomes by regression with random effects." \emph{Statistics in Medicine}. \bold{17}(22):2537-2550.

Ritz, J., E. Demidenkob, et al. (2008). "Multivariate meta-analysis for data consortia, individual patient meta-analysis, and pooling projects." \emph{Journal of Statistical Planning and Inference}. \bold{139}(7):1919-1933.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\seealso{ 
See \code{\link{qtest}} for the generic method function. See \code{\link{mvmeta-package}} and \code{\link{mvmeta}} for an overview of the package and modelling framework.
}

\examples{
# RUN THE MODEL AND SUMMARIZE THE RESULTS
model <- mvmeta(y=berkey98[4:5],S=berkey98[6:8])

# MULTIVARIATE COCHRAN Q TEST FOR HETEROGENEITY
test <- qtest(model)
print(test,digits=2)

}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
\keyword{htest}
