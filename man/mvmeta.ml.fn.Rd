\name{mvmeta.ml.fn}
\alias{mvmeta.ml.fn}
\alias{mvmeta.reml.fn}
\alias{mvmeta.ml.gr}
\alias{mvmeta.reml.gr}
\alias{mvmeta.igls}

\title{ Likelihood functions for mvmeta Models }

\description{
These functions computes the value of the log-likelihood and the related vectors of first partial derivatives for random-effects multivariate and univariate meta-analysis and meta-regression, in terms of model parameters. They are meant to be used internally and not directly run by the users.
}

\usage{

mvmeta.ml.fn(par, Xlist, ylist, Slist, nalist, k, m, nall)
mvmeta.ml.gr(par, Xlist, ylist, Slist, nalist, k, m, nall)

mvmeta.reml.fn(par, Xlist, ylist, Slist, nalist, k, m, nall)
mvmeta.reml.gr(par, Xlist, ylist, Slist, nalist, k, m, nall)

mvmeta.igls(Psi, Xlist, ylist, Slist, nalist, k, m)
}

\arguments{
Assuming a meta-analysis or meta-regression based on \eqn{m} studies, \eqn{k} outcomes and \eqn{p} predictors:

  \item{par }{ a \eqn{k(k+1)/2}-dimensional vector representing the upper triangle entries of the Cholesky decomposition of \code{Psi}, taken by row.}
  \item{Psi }{ a \eqn{k \times k}{k x k} matrix representing the current estimate of the between-study (co)variance matrix.}
  \item{Xlist }{ a \eqn{m}-dimensional list of study-specific design matrices for the fixed-effects part of the model. Rows corresponding to missing outcomes have been excluded.}
  \item{ylist }{ a \eqn{m}-dimensional list of study-specific of vectors of estimated outcomes. Entries corresponding to missing outcomes have been excluded.}
  \item{Slist }{ a \eqn{m}-dimensional list of within-study (co)variance matrices of estimated outcomes. Rows and columns corresponding to missing outcomes have been excluded.}
  \item{nalist }{ a \eqn{m}-dimensional list of \eqn{k}-dimensional study-specific logical vectors, identifying missing outcomes.}
  \item{k, m, nall }{ numeric scalars: number of outcomes, number of studies included in estimation (equal to the length of lists above), number of observations (excluding missing).}

}

\details{
These functions are called internally by the fitting functions \code{\link{mvmeta.ml}} and \code{\link{mvmeta.reml}} to perform iterative optimization algorithms for estimating random effects meta-analytical models.

The maximization of the (restricted) likelihood starts with few runs of an iterative generalized least square algorithm implemented in \code{mvmeta.igls}, which can be regarded as a fast and stable way to get starting values close to the maximum for the Quasi-Newton iterative algorithm, implemented in \code{\link{optim}}. This function requires the functions to compute the value of the (restricted) likelihood in terms of the model parameters, and the vector of first partial derivatives, provided by the functions \code{mvmeta.ml.fn}-\code{mvmeta.reml.fn} and \code{mvmeta.ml.gr}-\code{mvmeta.reml.gr}, respectively.

These functions actually specify the \emph{profiled} version of the (restricted) likelihood, basically expressed only in terms of the components of \code{Psi}, while the estimate of the fixed-effects coefficients is provided at each iteration by the internal function \code{.gls}, based on the current value of \code{Psi}. At convergence, the value of this profiled version is identical to the full (restricted) likelihood. This approach is computationally efficient, as it reduces the number of parameters in the optimization routine, especially for meta-regression models.

The parameters of \code{Psi} are expressed in terms of its Cholesky decomposition, in order to assure the positive-definiteness of the estimated between-study (co)variance matrix. A second Cholesky decomposition is then performed on the marginal (co)variance matrix in order to re-express the problem as standard least square equations, an approach which speeds up the computation of matrix inverses and determinants. These equations are finally solved through a QR decomposition, which guarantees stability. More details are provided in the references below.
}

\value{
\code{mvmeta.ml.fn} and \code{mvmeta.reml.fn} return the value of the (restricted) log-likelihood for a given set of parameters in \code{par}. \code{mvmeta.ml.gr} and \code{mvmeta.reml.gr} return instead the related vector of first partial derivatives. \code{mvmeta.igls} returns an updated estimate of \code{Psi} given its initial value or the value at previous iteration.
}

\references{
Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. Epub ahead of print (doi 10.1002/sim.5471). [Freely available \href{http://onlinelibrary.wiley.com/doi/10.1002/sim.5471/pdf}{\bold{here}}].

Goldstein H (1986). Multilevel mixed linear model analysis using iterative generalized least squares. \emph{Biometrika}. \bold{73}(1):43.

Lindstrom MJ and Bates DM (1988). Newton-Raphson and EM algorithms for linear mixed-effects models for repeated-measures data. \emph{Journal of the American Statistical Association}. \bold{83}(404):1014--1022.

Pinheiro JC and Bates DM (2000). Mixed-Effects Models in S and S-PLUS. New York, Springer Verlag.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
As stated earlier, these functions are called internally by \code{\link{mvmeta.ml}} and \code{\link{mvmeta.reml}}, and are not meant to be used directly. In particular, their code does not contain any check on the arguments provided, which are expected in specific formats. They are however exported in the namespace and documented for completeness.
}

\seealso{
See \code{\link{mvmeta.fit}} and the \code{\link[=mvmeta.ml]{fitting functions}} for additional info on the fitting procedures. See \code{\link{chol}} and \code{\link{qr}} for info on the Cholesky and QR decomposition.
See \code{\link{mvmeta-package}} for an overview of the package and modelling framework.
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
