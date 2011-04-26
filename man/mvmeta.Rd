\name{mvmeta}
\alias{mvmeta}


\title{ Fitting Multivariate Meta-Analysis and Meta-Regression }

\description{
This function fits fixed and random-effects multivariate meta-analysis and meta-regression models. See \code{\link{mvmeta-package}} for an overview of this modelling framework.
}

\usage{
mvmeta(y, S, X=NULL, method="reml", lab, cen=FALSE, na.action, ...)
}

\arguments{
In multivariate meta-analysis, \eqn{k} outcome parameters \eqn{y}, estimated in \eqn{m} studies, are pooled accounting for their (co)variance matrix \eqn{S}. Optionally, \eqn{p-1} predictors (plus intercept) may be included  in order specify a multivariate meta-regression model. The arguments of the function are:

  \item{y }{ series of \eqn{k} estimated outcome parameters for each one of the \eqn{m} studies. Accepted formats are a \eqn{m} x \eqn{k} matrix or dataframe; or alternatively a \eqn{m}-length list of \eqn{k}-length vectors.}
  \item{S }{ series of within-study (co)variance matrices of the estimated outcome parameters for each one of the \eqn{m} studies. Accepted formats are a \eqn{m}-length list of \eqn{k} x \eqn{k} matrices; a tri-dimensional \eqn{k} x \eqn{k} x \eqn{m} array; or a \eqn{m} x \eqn{[k(k+1)/2]} matrix or dataframe where each row represents the entries of the lower triangle of the related (co)variance matrix, taken by column.}
  \item{X }{ series of \eqn{p-1} study-level predictors (excluding intercept) for each one of the \eqn{m} studies. Accepted formats are a \eqn{m} x \eqn{(p-1)} matrix or dataframe; or alternatively a \eqn{m}-length list of \eqn{(p-1)}-length vectors. For a single predictor, aslo a vector is accepted.}
  \item{method }{ estimation method: \code{"fixed"} for fixed-effects models, \code{"ml"} or \code{"reml"} for random-effects models through (restricted) maximum likelihood.}
  \item{lab }{ list of label vectors: \code{mlab} for the \eqn{m} studies, \code{klab} for the \eqn{k} outcome parameters, \code{plab} for the \eqn{(p-1)} predictors. Set to default values if not provided.}
  \item{cen }{ centering values for the predictors in meta-regression. Accepted formats are a \eqn{(p-1)}-length vector, or logical. If \code{TRUE}, centered on predictors means, if \code{FALSE} (default) left uncentered.}
  \item{na.action }{ a function which indicates what should happen when the data contain \code{NA}s. The default is set by the \code{na.action} setting of \code{\link{options}}, usually \code{\link{na.omit}}. Value \code{\link{na.exclude}} can be useful.}
  \item{\dots }{ further argument passed to functions called internally. See Details.}

}

\details{
Fixed-effects models are estimated through generalized least square given the within-study (co)variance matrices \eqn{S}. Random-effects models are estimated through profiled (restricted) maximum likelihood, with quasi-Newton iterative algorithms. In this approach the (restricted) likelihood is maximized exclusively with respect to parameters of the between-study (co)variance matrix \eqn{Psi}, with coefficients \eqn{beta} re-computed at each iteration. Additional information on the estimation procedures are given in pages of \code{\link[mvmeta:mvmeta.ml]{estimation algorithms}}.

Missing valus are allowed in \code{y}, \code{S} and \code{X}. In the latter case, the related study is entirely excluded from estimation. Instead, if \code{y} and \code{S} are partially provided, the function excludes in each study the outcome parameters for which the estimate or one of the (co)variance terms are missing.

In meta-regression models, categorical predictors (factors with dummy parameterization) are not correctly handled if included in the argument \code{X}. Eventually, they must be manually expanded in the correct matrix form, for example using \code{\link{model.matrix}}. In this case, the intercept term must be excluded. This limitation is likely to be fixed in next releases.

Centering values specified by \code{cen} may be useful for interpretation purposes in meta-regression models. Specifically, the intercept coefficients are interpreted as the estimate of the outcome parameters for values of the predictors equal to the centering points.

Additional arguments specified by \dots are passed to the function \code{\link{optim}}, called internally by \code{\link{mvmeta}} for random-effects models, and may be used to control the optimization procedure. See the related help page for a list of specific arguments.
}

\value{
A list object of class \code{"mvmeta"} with the following components:
  \item{y }{ A list of \eqn{k}-length vectors of estimated outcome parameters rearranged from the argument above, with labels.}
  \item{S }{ A list of \eqn{k} x \eqn{k} (co)variance matrices rearranged from the argument above, with labels.}
  \item{X }{ A \eqn{m} x \eqn{(p-1)} predictors matrix rearranged from the argument above, with labels.}
  \item{beta, vcov }{ estimated coefficients and related (co)variance matrix for the fixed part of the model.}
  \item{Psi }{ estimated between-study (co)variance matrix.}
  \item{method }{ as above.}
  \item{logLik }{ (restricted) maximum likelihood.}
  \item{df }{ list with the following scalar components: \code{nobs} (number of observations, excluding missing values), \code{fixed} (number of estimated coefficients), \code{random} (number of estimated (co)variance parameters).}
  \item{dim }{ list with the following scalar components: \code{m} (number of studies included in estimation, which could be lower than the total number in the presence of missing values), \code{k} (number of parameters), \code{p} (number of coefficients for each outcome parameter).}
  \item{lab }{ list of label vectors as above.}
  \item{cen }{ as above.}
  \item{na.action }{ as above.}
  \item{convergence }{ An integer code. 0 indicates successful completion, 1 that the iteration limit has been reached without convergence.}
The method functions \code{\link{print.mvmeta}} and \code{\link{summary.mvmeta}} return a summary of the results.
}

\references{
Jackson, D., Riley, R., White, I. R. Multivariate meta-analysis: Potential and promise. \emph{Statistics in Medicine}. Online First. DOI: 10.1002/sim.4172.

Berkey, C. S., D. C. Hoaglin, et al. (1998). "Meta-analysis of multiple outcomes by regression with random effects." \emph{Statistics in Medicine}. \bold{17}(22):2537-2550.

White, I. R. (2009). "Multivariate random-effects meta-analysis." \emph{Stata Journal}. \bold{9}(1):40-56.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
The current implementation of the function does not allow any pre-specified structure on the matrix \eqn{Psi}, and requires the estimation of  \eqn{k(k+1)/2} (co)variance parameters. Simpler structure will be included in the next releases. In addition, meta-regression models specify \eqn{p} coefficients for each one of the \eqn{k} outcome parameters. Future improvements will offer the option to define different equation formulae for each outcome, for example excluding the intercept or associating to each outcome a subset of the other \eqn{(p-1)} predictors.
}

\section{Warnings}{
This release of the package \pkg{mvmeta} has been tested with different simulated and real datasets. The functions generally perform very well under several scenarios, in comparison with alternative software implementations. However, bugs and bad performances under un-tested conditions may not be excluded. Please report any error or unexpected behaviour to the e-mail address below.
}

\seealso{ 
\code{\link{mvmeta-package}} for an overview of this modelling framework. Pages of \code{\link[mvmeta:mvmeta.ml]{estimation algorithms}} for additional information on the estimation procedures.
}

\examples{
# INSPECT THE DATA
berkey98

### FIRST EXAMPLE
### MULTIVARIATE META-ANALYSIS
### ESTIMATED THROUGH REML

# RUN THE MODEL AND SUMMARIZE THE RESULTS
model <- mvmeta(y=berkey98[4:5],S=berkey98[6:8])
summary(model)
# AIC VALUE
AIC(model)
# BLUP ESTIMATES AND 90% PREDICTION INTERVALS, AGGREGATED BY OUTCOME
blup(model,pi=TRUE,aggregate="y",pi.level=0.90)

### SECOND EXAMPLE
### MULTIVARIATE META-REGRESSION
### ESTIMATED THROUGH ML, CENTERED ON PREDICTOR MEAN

# RUN THE MODEL AND SUMMARIZE THE RESULTS
model <- mvmeta(y=berkey98[4:5],S=berkey98[6:8],
	X=berkey98$pubyear,cen=TRUE,method="ml")
summary(model)
# MULTIVARIATE COCHRAN Q TEST FOR RESIDUAL HETEROGENEITY
qtest(model)
# PREDICTED AVERAGED OUTOCOMES AND STANDARD ERRORS FROM YEAR 1985 TO 1989
predict(model,newdata=1985:1989,se=TRUE)
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}

