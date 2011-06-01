\name{mvmeta}
\alias{mvmeta}


\title{ Fitting Multivariate Meta-Analysis and Meta-Regression }

\description{
This function fits fixed and random-effects multivariate meta-analysis and meta-regression models. See \code{\link{mvmeta-package}} for an overview of this modelling framework.
}

\usage{
mvmeta(formula, S, data, method="reml", lab, contrasts, na.action, control)
}

\arguments{
In multivariate meta-analysis, \eqn{k} outcome parameters \eqn{y}, estimated in \eqn{m} studies, are pooled accounting for their (co)variance matrix \eqn{S}. Optionally, \eqn{p} predictors (usually with intercept) may be included  in order specify a multivariate meta-regression model. The arguments of the function are:

  \item{formula }{ an object of class \code{"\link{formula}"} (or one that can be coerced to that class) offering a symbolic description of the model to be fitted. Alternatively, for meta-analysis with no predictor, a single matrix -type object. Terms must be matrix or vectors objects, optionally provided in the \code{data} argument below. See Details.}
  \item{S }{ series of within-study (co)variance matrices of the estimated outcome parameters for each one of the \eqn{m} studies. Accepted formats are a \eqn{m}-length list of \eqn{k} x \eqn{k} matrices; a tri-dimensional \eqn{k} x \eqn{k} x \eqn{m} array; or a \eqn{m} x \eqn{[k(k+1)/2]} matrix or dataframe where each row represents the entries of the lower triangle of the related (co)variance matrix, taken by column. Optionally, terms may be provided in the \code{data} argument below.}
  \item{data }{ an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in \code{formula}. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{mvmeta} is called.}
  \item{method }{ estimation method: \code{"fixed"} for fixed-effects models, \code{"ml"} or \code{"reml"} for random-effects models through (restricted) maximum likelihood.}
  \item{lab }{ list of label vectors: \code{mlab} for the \eqn{m} studies, \code{klab} for the \eqn{k} outcome parameters, \code{plab} for the \eqn{(p)} predictors (including intercept and allowing for factor expansion). Optionally, \code{mlab} may be provided in the \code{data} argument above. Set to default values if not provided.}
  \item{contrasts }{ an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.}
  \item{na.action }{ a function which indicates what should happen when the data contain \code{NA}s. The default is set by the \code{na.action} setting of \code{\link{options}}, usually \code{\link{na.omit}}. Value \code{\link{na.exclude}} can be useful. See Details.}
  \item{control }{ a list of further argument passed to functions called internally. See Details.}

}

\details{
Multivariate meta-analysis is specified with the formula \code{y ~ 1}, where the left-hand side is a matrix (or vector for univariate analysis), optionally of form \code{cbind(y1,...,yk)} with terms stored in \code{data}. Alternatively, matrix or vector objects are allowed, and the formula is retrieved internally adding \code{~ 1}. In (multivariate) meta-regression, other terms are added in the right-hand side of the formula. Factors, variable transformations and interactions are allowed, following the standard formula specification. See \code{\link{formula}} for further details. 

Fixed-effects models are estimated through generalized least square given the within-study (co)variance matrices \eqn{S}. Random-effects models are estimated through profiled (restricted) maximum likelihood, with quasi-Newton iterative algorithms. In this approach the (restricted) likelihood is maximized exclusively with respect to parameters of the between-study (co)variance matrix \eqn{Psi}, with coefficients \eqn{beta} re-computed at each iteration. Additional information on the estimation procedures are given in pages of \code{\link[mvmeta:mvmeta.ml]{estimation algorithms}}.

Missing valus are allowed in both sides of \code{formula} and in \code{S}. In the case of missing predictors (right-hand side of \code{formula}), the related study is entirely excluded from estimation. Instead, if outcomes or their (co)variance matrix in \code{S} are partially provided, the function excludes in each study the outcome parameters for which the estimate or one of the (co)variance terms are missing.

Additional arguments specified by \code{control} are passed to the function \code{\link{optim}}, called internally by \code{mvmeta} for random-effects models, and may be used to control or inspect the optimization procedure. See the related help page for a list of specific arguments.
}

\value{
A list object of class \code{"mvmeta"} with the following components:

  \item{beta, vcov }{ estimated coefficients and related (co)variance matrix for the fixed part of the model.}
  \item{Psi }{ estimated between-study (co)variance matrix.}
  \item{logLik }{ (restricted) maximum likelihood.}
  \item{df }{ list with the following scalar components: \code{nobs} (number of observations, excluding missing values), \code{fixed} (number of estimated coefficients), \code{random} (number of estimated (co)variance parameters).}
  \item{dim }{ list with the following scalar components: \code{m} (number of studies included in estimation, which could be lower than the total number in the presence of missing values), \code{k} (number of parameters), \code{p} (number of coefficients for each outcome parameter).}
  \item{lab }{ list of label vectors as above.}
  \item{int }{ logical switch indicating if an intercept is included in the model formula.}
  \item{call }{ the function call.}
  \item{model }{ the data frame used for estimation including outcomes and predictors variables, with attributes.}
  \item{convergence }{ An integer code. 0 indicates successful completion, 1 that the iteration limit has been reached without convergence.}

Other arguments are returned as specified in the Usage section. The method functions \code{\link{print.mvmeta}} and \code{\link{summary.mvmeta}} return a summary of the results.
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
model <- mvmeta(cbind(PD,AL),S=berkey98[6:8],data=berkey98,lab=list(mlab=trial))
summary(model)
# AIC VALUE
AIC(model)
# BLUP ESTIMATES AND 90% PREDICTION INTERVALS, AGGREGATED BY OUTCOME
blup(model,pi=TRUE,aggregate="y",pi.level=0.90)

### SECOND EXAMPLE
### MULTIVARIATE META-REGRESSION
### ESTIMATED THROUGH ML

# RUN THE MODEL AND SUMMARIZE THE RESULTS
model <- mvmeta(cbind(PD,AL)~pubyear,S=berkey98[6:8],data=berkey98,
	lab=list(mlab=trial),method="ml")
summary(model)
# MULTIVARIATE COCHRAN Q TEST FOR RESIDUAL HETEROGENEITY
qtest(model)
# PREDICTED AVERAGED OUTOCOMES AND STANDARD ERRORS FROM YEAR 1985 TO 1989
newdata <- data.frame(pubyear=1985:1989)
predict(model,newdata,se=TRUE)
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}

