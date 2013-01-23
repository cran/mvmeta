\name{mvmeta-package}
\alias{mvmeta-package}
\docType{package}

\title{
Multivariate Meta-Analysis and Meta-Regression
}

\description{
The package \pkg{mvmeta} consists of a collection of functions to perform fixed and random-effects multivariate and univariate meta-analysis and meta-regression in \R.
}

\section{Modelling framework}{
Multivariate meta-analytical models represent an extension of the standard univariate techniques, where estimates of a single effect size, here defined generally as outcome, are pooled across studies. In multivariate meta-analysis, estimates of multiple outcomes are combined while accounting for their correlation. Multivariate meta-regression also models such multivariate distribution in terms of study-level predictors. These statistical tools have been originally proposed to model multiple endpoints in clinical trials. Applications and methodological developments are currently proposed also for network meta-analysis (indirect treatment comparison) and for the meta-analysis of multi-parameter associations.

Similarly to univariate methods, fixed-effects models do not assume heterogeneity among studies, and the estimates are conditional on the set of studies collected in the meta-analysis, based on a multivariate weighted average of study-level estimates. Random-effects meta-analysis, instead, allows a degree of heterogeneity among studies, assuming the (true but unobserved) study-specific outcomes as randomly sampled from a multivariate normal distribution of studies. Inference from these models may therefore be extended to other unmeasured studies assumed to belong to the same (usually hypothetical) population.

Estimation and interpretation exploit here the framework of linear mixed models. The fixed part of the model provides an estimate of fixed effects, which represent the population-averaged outcomes and, in the case of multivariate meta-regression, are defined by a set of coefficients associated to study-level predictors. The random part of the model describes the deviation from the population averages, estimating the components of a between-study (co)variance matrix. Assuming \eqn{k} outcomes \eqn{\boldsymbol{y}_i}{yi} estimated in each of \eqn{i = 1, \ldots, m} studies, and related to \eqn{p} study-level predictors \eqn{\boldsymbol{x}_i}{xi}, random-effects multivariate meta-regression models can be generally described with:

\deqn{\boldsymbol{y}_i \sim \mbox{N}_k(\boldsymbol{X}_i \boldsymbol{\beta}, \boldsymbol{S}_i + \Psi)}{yi ~ Nk(Xi\beta , Si + \Psi)}

Here the outcomes \eqn{\boldsymbol{y}_i}{yi} are assumed to be sampled from a multivariate normal distribution of order \eqn{k}. Their distribution is centred on \eqn{\boldsymbol{X}_i \boldsymbol{\beta}}{Xi\beta}, with \eqn{\boldsymbol{X}_i}{Xi} as a \eqn{k \times kp}{k x kp} design matrix and \eqn{\boldsymbol{\beta}}{\beta} the vector of fixed-effects coefficients. The marginal \eqn{k \times k}{k x k} (co)variance matrix \eqn{\boldsymbol{\Sigma}_i = \boldsymbol{S}_i + \boldsymbol{\Psi}}{\Sigmai = Si + \Psi} is given by the sum of within (assumed known) and between-study (co)variance matrices \eqn{\boldsymbol{S}_i}{Si} and \eqn{\boldsymbol{\Psi}}{\Psi}, respectively.

Other models are taken as special cases of that above. In multivariate meta-analysis, \eqn{\boldsymbol{X}_i}{Xi} becomes an identity matrix with \eqn{p=1}, and \eqn{\boldsymbol{\beta}}{\beta} reduces to \eqn{k} intercepts, interpreted as the population-averaged outcomes. For \eqn{k=1}, the model reduces to the standard univariate meta-analysis or meta-regression. In fixed-effects meta-analytic models, \eqn{\boldsymbol{\Psi}}{\Psi} is assumed not to exist, and the variability between studies is due exclusively to the within-study estimation error.
}

\section{Estimation methods}{
The aim is to estimate the coefficients \eqn{\boldsymbol{\beta}}{\beta} and, for random-effects models, the components of the between-study (co)variance matrix \eqn{\boldsymbol{\Psi}}{\Psi}. If the same linear predictor of \eqn{p} terms \eqn{\boldsymbol{x}_i}{xi} measured in each study is specified for all the outcomes (the only option available at the moment), the dimension of \eqn{\boldsymbol{\beta}}{\beta} is \eqn{kp}. The parameters for the random part depend on the chosen structure of the between-study (co)variance matrix \eqn{\boldsymbol{\Psi}}{\Psi}, with \eqn{k(k+1)/2} parameters for an unstructured form (the only option available at the moment). 

The estimation options available in the \pkg{mvmeta} package are:

\itemize{
  \item \bold{Fixed-effects}
  \item \bold{Maximum likelihood} (ML)
  \item \bold{Restricted maximum likelihood} (REML)
  \item \bold{Method of moments}
}

The fixed-effects model is fitted through generalized least squares (GLS), assuming the (co)variance structure, composed by the within-study error only, as completely known. Likelihood-based methods exploit on a quasi-Newton iterative algorithm, profiled on the parameters for the (co)variance components, where the estimates for the fixed effects are estimated through GLS at each iteration. The method of moments is a multivariate extension of the traditional estimator used in univariate models.

Likelihood-based random-effects models are computationally intensive, and convergence can be slow for high-dimensional models (with a high number of outcomes). However, fit criteria and inferential test derived from likelihood theory, such as AIC and likelihood ratio test, can be derived. In contrast, the semi-parametric approach based on method of moments is non-iterative (meaning sometimes faster), and do not require the usual normality assumption for the distribution of the random effects, although statistics for model comparison are not available. Both approaches apply constraints for the between-study (co)variance matrix \eqn{\boldsymbol{\Psi}}{\Psi} to be positive semi-definite. Further details on estimation methods are given in the related help pages.
}

\section{Functions and data included in the package}{
The main function in the package is \code{\link{mvmeta}}, which performs the various models illustrated above. This function resembles standard regression functions in \R, and specifies the model through a regression formula. The function returns a list object of class \code{"mvmeta"} (see \code{\link{mvmetaObject}}).

The estimation is carried out through \code{\link{mvmeta.fit}}, a wrapper which prepares the data and calls internal functions for fitting the models. Specifically, \code{\link{mvmeta.fixed}} is applied for fixed-effects models, \code{\link{mvmeta.ml}} and \code{\link{mvmeta.reml}} are called for estimating random-effects models through (restricted) maximum likelihood, while \code{\link{mvmeta.mm}} is used for random-effects models estimated through the method of moments. For likelihood-based methods, iterative optimization algorithms for maximizing the (restricted) likelihood exploit the functions \code{\link{mvmeta.ml.fn}}-\code{\link{mvmeta.ml.gr}} and \code{\link{mvmeta.reml.fn}}-\code{\link{mvmeta.reml.gr}}. Starting values are computed with the function \code{\link{mvmeta.igls}}. Fitting parameter options are set by \code{\link{mvmeta.control}}.

Method functions are available for objects of class \code{"mvmeta"} (see \code{\link{mvmetaObject}} for a complete list). \code{\link[=summary.mvmeta]{summary}} produces a list of class \code{"summary.mvmeta"} for summarizing the fit of the model and providing additional results. The method function \code{\link[=predict.mvmeta]{predict}} computes predicted values, optionally for a set of new values of the predictors. \code{\link[=blup.mvmeta]{blup}} gives the (empirical) best linear unbiased predictions for the set of studies used for estimation. Other default or specific method functions for regression can be used on objects of class \code{"mvmeta"}, such as \code{\link{fitted}} and \code{\link{residuals}}, \code{\link[=logLik.mvmeta]{logLik}}, \code{\link{AIC}} and \code{\link{BIC}}, among others.

Methods for \code{\link[=model.frame.mvmeta]{model.frame}} and \code{\link[=model.matrix.mvmeta]{model.matrix}} are used to extract and construct the model frame and the design matrix of the regression meta-analytical model, respectively. Methods for \code{\link[=na.omit.data.frame.mvmeta]{na.omit}} and \code{\link[=na.exclude.data.frame.mvmeta]{na.exclude}} help handle correctly missing values.

Simulations can be produced using the function \code{\link{mvmetaSim}} and the method function \code{\link[=simulate.mvmeta]{simulate}}, which return one or multiple sets of simulated outcomes for a group of studies.

The method function \code{\link{qtest.mvmeta}} (producing an object with class of the same name) performs the (multivariate) Cochran Q test for (residual) heterogeneity, both on the overall multivariate distribution and on each single outcome. The generic method function is \code{\link{qtest}}.

Printing functions for the objects of classes defined above are also provided. Other internal functions of the package (with name starting with a dot) are not exported in the namespace. For users interested in getting into details of the package structure, the function \code{\link{getAnywhere}} can be of help to access the code of the internal functions.

The dataset \code{\link{berkey98}} includes the results for 2 outcomes in 5 trials on periodontal diseases, and is used in the examples.
}

\section{Future developments}{
The package \pkg{mvmeta} will hopefully experience substantial changes and improvements in the next releases. Functions to compute residuals and other model checking methods will be included. If possible, more parsimonious structures (e.g. diagonal, compound-symmetry) for the between-study covariance matrix will be allowed, reducing the number of parameters in high-dimensional models.

A list of changes included in the current and previous versions can be found by typing:

\code{file.show(system.file("ChangeLog",package="mvmeta"))}

General information on the development and applications of the \pkg{mvmeta} package and on the modelling framework of multivariate meta-analysis, together with an updated version of the \R{} scripts for running the examples in published papers, can be found at \bold{\href{http://www.ag-myresearch.com}{www.ag-myresearch.com}}.
}

\note{
The function included in the package \pkg{mvmeta} may be applied also to perform standard univariate meta-analysis and meta-regression. However, alternative packages provides a more exhaustive and efficient set of functions. See for example the packages \pkg{metafor}, \pkg{meta} and \pkg{rmeta}, among others.

The package \pkg{mvtmeta} performs fixed-effects and random-effects multivariate meta-analysis using the same method of moments estimator adopted here, although without allowing for missing outcomes or meta-regression.

Use \code{citation("mvmeta")} to cite this package.
}

\section{Warnings}{
This release of the package \pkg{mvmeta} has been tested with different simulated and real datasets. The functions generally perform well under several scenarios, and comparisons with alternative software implementations show good agreement. However, bugs and bad performances under un-tested conditions may not be excluded. Please report any error or unexpected behaviour to the e-mail address below.
}

\author{
Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}
}

\references{
Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/statmed2012.html}{\bold{here}}].

Jackson D, Riley R, White IR (2011). Multivariate meta-analysis: Potential and promise. \emph{Statistics in Medicine}. \bold{30}(20);2481--2498.

Berkey, CS, Anderson JJ, Hoaglin DC (1996). Multiple-outcome meta-analysis of clinical trials. \emph{Statistics in Medicine}. \bold{15}(5):537--547.

Berkey, CS, Hoaglin DC, et al. (1998). Meta-analysis of multiple outcomes by regression with random effects. \emph{Statistics in Medicine}. \bold{17}(22):2537--2550.

White IR (2009). Multivariate random-effects meta-analysis. \emph{Stata Journal}. \bold{9}(1):40--56.

White IR (2011). Multivariate random-effects meta-regression: updates to mvmeta. \emph{Stata Journal}. \bold{11}(2):255--270.

Jackson D, White IR, Riley RD (2013). A matrix based method of moments for fitting the multivariate random effects model for meta-analysis and meta-regression. \emph{Biometrical Journal}. To be published.

Chen H, Manning AK, Dupuis J (2012). A method of moments estimator for random effect multivariate meta-analysis. \emph{Biometrics}. \bold{68}(4):1278-1284.

van Houwelingen HC, Arends LR, et al. (2002). Advanced methods in meta-analysis: multivariate approach and meta-regression. \emph{Statistics in Medicine}. \bold{21}(4):589--624.

Nam IS, Mengersen K, et al. (2003). Multivariate meta-analysis. \emph{Statistics in Medicine}. \bold{22}(14):2309--2333.

Arends LR, Voko Z, Stijnen T (2003). Combining multiple outcome measures in a meta-analysis: an application. \emph{Statistics in Medicine}. \bold{22}(8):1335--1353.

Ritz J, Demidenkob E, Spiegelman G (2008). Multivariate meta-analysis for data consortia, individual patient meta-analysis, and pooling projects. \emph{Journal of Statistical Planning and Inference}. \bold{139}(7):1919--1933.
}

\keyword{package}
