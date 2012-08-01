\name{mvmeta}
\alias{mvmeta}
\alias{mvmeta.fit}


\title{ Fitting Multivariate Meta-Analysis and Meta-Regression Models }

\description{
The function \code{mvmeta} performs fixed and random-effects multivariate and univariate meta-analysis and meta-regression, with various estimation methods. The function \code{mvmeta.fit} is a wrapper for actual fitting functions based on different estimation methods, usually called internally. See \code{\link{mvmeta-package}} for an overview.
}

\usage{
mvmeta(formula, S, data, subset, method="reml", model=TRUE, contrasts=NULL,
    offset, na.action, control=list())
    
mvmeta.fit(X, y, S, offset=NULL, method="reml", control=list())
}

\arguments{
Assuming a meta-analysis or meta-regression based on \eqn{m} studies, \eqn{k} outcomes and \eqn{p} predictors:

  \item{formula }{ an object of class \code{"\link{formula}"} (or one that can be coerced to that class) offering a symbolic description of the linear predictor of the model to be fitted to each outcome. Alternatively, for meta-analysis with no predictor, a single  vector (for univariate models) or matrix-type object (for multivariate models). Terms in \code{formula} must be vector or matrix-type objects, optionally provided in the \code{data} argument below. See Details.}
  \item{X }{ a \eqn{m \times p}{m x p} design matrix containing the study-specific predictors. Usually produced internally by \code{mvmeta} from \code{formula} above.}
  \item{y }{ a \eqn{m}-dimensional vector (for univariate models) or \eqn{m \times k}{m x k} matrix (for multivariate models) of outcomes. Usually produced internally by \code{mvmeta} from \code{formula} above.}
  \item{S }{ series of within-study (co)variance matrices of the estimated outcome for each one of the \eqn{m} studies. Accepted formats by \code{mvmeta} are a \eqn{m}-dimensional list of \eqn{k \times k}{k x k} matrices; a tri-dimensional \eqn{k \times k \times m}{k x k x m} array; or a \eqn{m \times k(k+1)/2}{m x k(k+1)/2} matrix or dataframe where each row represents the vectorized entries of the lower triangle of the related (co)variance matrix, taken by column (see \code{\link{xpndMat}}). Optionally, terms may be provided in the \code{data} argument below. \code{mvmeta.fit} accepts only the last option.}
  \item{data }{ an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in \code{formula}. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{mvmeta} is called.}
  \item{subset }{ an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{method }{ estimation method: \code{"fixed"} for fixed-effects models, \code{"ml"} or \code{"reml"} for random-effects models through (restricted) maximum likelihood. If \code{"model.frame"}, the model frame is returned. See \code{\link{lm}} or \code{\link{glm}}.}
  \item{model }{ a logical value indicating whether the model frame should be included as a component of the returned value. See the \code{\link[=model.frame.mvmeta]{model.frame}} method function.}
  \item{contrasts }{ an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix}}.}
  \item{offset }{ optionally, a \eqn{m}-dimensional numeric vector used to specify an a priori known component in the linear predictor. One or more \code{\link{offset}} terms can be included in the formula instead or as well. See \code{\link{model.offset}}.}
  \item{na.action }{ a function which indicates what should happen when the data contain \code{NA}s. Default to \code{na.action} setting of \code{\link{options}}, usually \code{\link[=na.omit.data.frame.mvmeta]{na.omit}}. \code{\link[=na.exclude.data.frame.mvmeta]{na.exclude}} can be useful. See details on \code{\link[=na.omit.data.frame.mvmeta]{missing values}} in \code{mvmeta}.}
  \item{control }{ list of parameters for controlling the fitting process. These are passed to  \code{\link{glm.control}} by \code{mvmeta.fit} to replace otherwise selected default values.}
}

\details{
The function \code{mvmeta} resembles standard regression functions in \R (see \code{\link{lm}} or \code{\link{glm}}). This function defines the design matrix and the vector (for univariate models) or matrix (for multivariate models) of outcomes, and calls the wrapper \code{mvmeta.fit} to perform the actual fitting. The latter prepares the data and calls specific \code{\link[=mvmeta.ml]{fitting functions}}, depending on the chosen \code{method}. Functions other than \code{mvmeta} are not expected to be called directly.

The model is specified through a regression formula. Simple meta-analysis is specified with the formula \code{y ~ 1}, where the left-hand side is a vector (in univariate models) or a matrix (in multivariate models), optionally of form \code{cbind(y1,...,yk)}, with terms stored in \code{data}. Alternatively, matrix or vector objects are allowed, and the formula is retrieved internally adding \code{~ 1}. In meta-regression, other terms are added in the right-hand side of the formula, defining the linear predictor common to all outcomes. Factors, variable transformations and interactions are allowed, following the standard formula specification. Labels are automatically retrieved from the objects in \code{formula}. See \code{\link{formula}} for further details. See \code{\link{lm}} or \code{\link{glm}} for info on the other arguments.

Fixed-effects models are estimated through generalized least square given the within-study (co)variance matrices \code{S}. Random-effects models are estimated through (restricted) maximum likelihood, with quasi-Newton iterative optimization algorithms provided by the function \code{\link{optim}} and based on specific \code{\link[=mvmeta.ml.fn]{likelihood functions}}.

Missing values are allowed in both sides of \code{formula}. In the case of missing predictors (right-hand side of \code{formula}), the related study is entirely excluded from estimation. In constrast, a study still contributes to estimation if at least outcome is non-missing. This behaviour is different than in standard regression functions such as \code{\link{lm}} or \code{\link{glm}}. Before the call to \code{mvmeta.fit}, studies matching such missing definition are removed from the the model frame. The missing pattern in \code{S} must be consistent with that in \code{y}. See further details on handling \code{\link[=na.omit.data.frame.mvmeta]{missing values}} in \code{mvmeta}.

The fitting procedure can be controlled through the additional terms specified in \code{control}, which are passed to the function \code{\link{mvmeta.control}}, called internally by \code{mvmeta.fit}. In the current versions, these are limited to the terms of the \code{control} argument of the function \code{\link{optim}} in random-effects models.
}

\value{
The \code{mvmeta} function typically returns a list object of class \code{"mvmeta"} representing the meta-analytical model fit, as described in \code{\link{mvmetaObject}}. When \code{method="data.frame"}, the model is not fitted and the model frame is returned, namely a data frame with special attributes (see the default method \code{\link{model.frame}}) and, in this case, the additional class \code{"data.frame.mvmeta"}.

The wrapper function \code{mvmeta.fit} is usually called internally in \code{mvmeta}, and returns an intermediate list object with some of the components expected in the \code{"mvmeta"} class.

Several method functions for regression objects are available, either default or specifically-written for the \code{"mvmeta"} class. See \code{\link{mvmetaObject}} for a complete list.
}

\references{
Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. Epub ahead of print (doi 10.1002/sim.5471). [Freely available \href{http://onlinelibrary.wiley.com/doi/10.1002/sim.5471/pdf}{\bold{here}}].

Jackson D, Riley R, White IR (2011). Multivariate meta-analysis: Potential and promise. \emph{Statistics in Medicine}. \bold{30}(20);2481--2498.

Berkey, CS, Hoaglin DC, et al. (1998). Meta-analysis of multiple outcomes by regression with random effects. \emph{Statistics in Medicine}. \bold{17}(22):2537--2550.

White IR (2009). Multivariate random-effects meta-analysis. \emph{Stata Journal}. \bold{9}(1):40--56.

White IR (2009). Multivariate random-effects meta-regression: updates to mvmeta. \emph{Stata Journal}. \bold{11}(2):255-270.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
In the current version, the same linear predictor must be set for all the outcomes, corresponding to that specified in \code{formula}. Also, only an unstructured form is supported for the between-study (co)variance matrix estimated as \code{Psi} (see \code{\link{mvmetaObject}}).
}

\seealso{ 
See additional info on the \code{\link[=mvmeta.ml]{fitting procedures}} and handling of \code{\link[=na.omit.data.frame.mvmeta]{missing values}} in \code{mvmeta}. See \code{\link{lm}} or \code{\link{glm}} for standard regression functions. See \code{\link{mvmeta-package}} for an overview of this modelling framework.
}

\examples{
# INSPECT THE DATA
berkey98

### FIRST EXAMPLE
### MULTIVARIATE META-ANALYSIS
### ESTIMATED THROUGH REML

# RUN THE MODEL
model <- mvmeta(cbind(PD,AL),S=berkey98[5:7],data=berkey98)

# SUMMARIZE THE RESULTS
summary(model)

# RESIDUALS AND FITTED VALUES
residuals(model)
fitted(model)

# LOG-LIKELIHOOD AND AIC VALUE
logLik(model)
AIC(model)

# BLUP ESTIMATES AND 90% PREDICTION INTERVALS, AGGREGATED BY OUTCOME
blup(model,pi=TRUE,aggregate="y",pi.level=0.90)


### SECOND EXAMPLE
### MULTIVARIATE META-REGRESSION
### ESTIMATED THROUGH ML

# RUN THE MODEL AND SUMMARIZE THE RESULTS
model <- mvmeta(cbind(PD,AL)~pubyear,S=berkey98[5:7],data=berkey98,method="ml")
summary(model)

# MULTIVARIATE COCHRAN Q TEST FOR RESIDUAL HETEROGENEITY
qtest(model)

# PREDICTED AVERAGED OUTOCOMES AND STANDARD ERRORS FROM YEAR 1985 TO 1989
newdata <- data.frame(pubyear=1985:1989)
predict(model,newdata,se=TRUE)

# MODEL FRAME AND MODEL MATRIX
model.frame(model)
model.matrix(model)


### THIRD EXAMPLE
### UNIVARIATE META-REGRESSION
### FIXED-EFFECTS MODEL

# RUN THE MODEL
model <- mvmeta(PD~pubyear,S=berkey98[,5],data=berkey98)

# SUMMARIZE THE RESULTS
summary(model)

# RESIDUALS AND FITTED VALUES
residuals(model)
fitted(model)

# COCHRAN Q TEST FOR RESIDUAL HETEROGENEITY
qtest(model)
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}

