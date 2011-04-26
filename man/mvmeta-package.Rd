\name{mvmeta-package}
\alias{mvmeta-package}
\docType{package}

\title{
Multivariate Meta-Analysis and Meta-Regression
}

\description{
The package \pkg{mvmeta} consists of a collection of functions to perform fixed and random-effects multivariate meta-analysis and meta-regression.
}

\section{Modelling framework}{
These models represent an extension of the standard meta-analityc techniques, where estimates of a single parameter (or outcome) are pooled across studies. In multivariate meta-analysis, multiple parameters are combined while accounting for their correlation. Multivariate meta-regression also models such multivariate distribution in terms of study-level predictors. These statistical tools have been proposed to model multiple endpoints in clinical trials, multi-parameter exposure-response relationships in observational studies, network meta-analysis, and may be potentially useful in other applications.

Similarly to univariate methods, fixed-effects models do not account for heterogeneity among studies, and the estimates are conditional on the set of studies collected in the meta-analysis. In practice, this approach provides a conditional inference on such set of studies, based on a multivariate weighted average of study-level estimates. Random-effects meta-analysis, instead, allows a degree of heterogeneity among studies, assuming the (true but unobserved) parameters as randomly sampled from a multivariate normal distribution of study population. Inference from these models may therefore be extended to other unmeasured studies belonging to the same (usually hypothetical) population.

Estimation and interpretation exploit the framework of linear mixed models. The fixed part of the model provides an estimate of fixed effects, representing the population-averaged outcome parameters and, in the case of multivariate meta-regression, the coefficients associated to study-level predictors. The random part of the model describes the deviation from the population averages, estimating the components of the between-study (co)variance matrix. These model may be summarized by the equation:

\deqn{y_i ~ N_k(kX_i beta, S_i + Psi)}

Here the \eqn{k} outcome parameters \eqn{y} estimated in each of the \eqn{i = 1, \ldots, m} studies are assumed to be sampled from a multivariate normal distribution of order \eqn{k}. The distribution is centered on \eqn{kX beta}, with \eqn{kX} a design matrix produced by the \eqn{p}-length vector of study-level predictors \eqn{x_i}, and the \eqn{kp}-length vector of coefficient \eqn{beta}. The covariance \eqn{Sigma = S + Psi} is given by the sum of within (assumed known) and between-study (co)variance matrices \eqn{S} and \eqn{Psi}, respectively. In multivariate meta-analysis with no predictors, \eqn{kX} is an identity matrix and \eqn{kX beta} reduces to \eqn{k} intercepts \eqn{beta}, representing here the vector of population-average outcome parameters. In fixed-effects meta-analytic models, \eqn{Psi} is assumed a zero matrix. The model aims to estimate the coefficient \eqn{beta} and, for random-effects models, the components of between-study (co)variance matrix \eqn{Psi}.

In multivariate meta-analysis, the \eqn{k} estimated coefficients \eqn{beta} are interpreted as the population-averaged outcome parameters, with \eqn{Psi} summarizing their (co)variance between study-specific realizations. In multivariate meta-regression, this variability (or part of it, for random-effects models) is explained by study-level predictors. In this case, the model specifies coefficients for the \eqn{p} predictors (including intercept) for each one of the \eqn{k} outcomes, with \eqn{kp} total fixed parameters in the model. 
}

\section{Functions included in the package}{
The main function in the package is \code{\link{mvmeta}}, which performs the various models illustrated above. The function accepts as arguments the estimates of parameters and their (co)variance matrix for \emph{m} studies, allowing different formats. Optionally, for multivariate meta-regression, study-level predictors may be included in the fixed part of the model. The function returns a list object of class \code{"mvmeta"}.

Estimation methods for random-effects models are based on (restricted) maximum likelihood, through quasi-Newton methods. The algorithms are provided by the internal functions \code{\link{mvmeta.ml}}, \code{\link{mvmeta.reml}}, \code{\link{mvmeta.ml.grad}}



The method function \code{\link[mvmeta:fitted.mvmeta]{fitted}} calculates fitted values for a given model, while \code{\link[mvmeta:predict.mvmeta]{predict}} compute predicted values for a set of new values of the predictors. \code{\link[mvmeta:blup.mvmeta]{blup}} gives the (shrunk) best linear unbiased prediction for the set of studies used for estimation.

The generic method function \code{\link{qtest}} and the method function \code{\link{qtest.mvmeta}} for \code{"mvmeta"} objects perform the (multivariate) Cochran Q test for (residual) heterogeneity.

The internal functions are not exported (not included in the package namespace), but nonetheless fully documented in the related help pages. Some other internal functions are worth being described: in particular, \code{\link{mvmeta.check}} performs coherence checks on the arguments of  \code{\link{mvmeta}}, while \code{\link{kXlistmk}} is called internally to build the design matrix for the fixed-effects part of the model.

The dataset \code{\link{berkey98}} includes the results for 2 outcomes in 5 trials on peridontal diseases, and is used in the examples.
}

\section{Future developments}{
The package \pkg{mvmeta} will hopefully experience substantial changes and improvements in the next releases. In particular, alternative estimation methods, such as multivariate extension of method of moments and Bayesian approaches, will be assessed and implemented. Functions to compute residuals and other model checking methods will be included. If possible, more parsimonious structures (e.g. diagonal, compound-symmetry) for the between-study covariance matrix will be allowed, reducing the number of parameters in high-dimensional models. In addition, different equations for different outcomes may be specified, assuming that some predictor is associated only with a subset of outcomes.

A list of changes included in the current and previous versions can be found by typing:

\code{file.show(system.file("ChangeLog",package="mvmeta"))}
}

\note{
The function included in the package \pkg{mvmeta} may be applied also to perform standard univariate meta-analysis and meta-regression. However, alternative packages provides a more exhaustive and efficient set of functions. See for example the packages \pkg{metafor}, \pkg{meta} and \pkg{rmeta}, among others.

Use \code{citation("mvmeta")} to cite this package.
}

\author{
Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}
}

\references{
Jackson, D., Riley, R., White, I. R. Multivariate meta-analysis: Potential and promise. \emph{Statistics in Medicine}. Online First. DOI: 10.1002/sim.4172.

Berkey, C. S., D. C. Hoaglin, et al. (1998). "Meta-analysis of multiple outcomes by regression with random effects." \emph{Statistics in Medicine}. \bold{17}(22):2537-2550.

White, I. R. (2009). "Multivariate random-effects meta-analysis." \emph{Stata Journal}. \bold{9}(1):40-56.

van Houwelingen, H. C., L. R. Arends, et al. (2002). "Advanced methods in meta-analysis: multivariate approach and meta-regression." \emph{Statistics in Medicine}. \bold{21}(4):589-624.

Nam, I. S., K. Mengersen, et al. (2003). "Multivariate meta-analysis." \emph{Statistics in Medicine}. \bold{22}(14):2309-2333.

Arends, L. R., Vok\'o, Z., et al. (2003). "Combining multiple outcome measures in a meta-analysis: an application." \emph{Statistics in Medicine}. \bold{22}(8):1335-1353.

Ritz, J., E. Demidenkob, et al. (2008). "Multivariate meta-analysis for data consortia, individual patient meta-analysis, and pooling projects." \emph{Journal of Statistical Planning and Inference}. \bold{139}(7):1919-1933.
}

\keyword{package}