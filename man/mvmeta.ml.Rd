\name{mvmeta.ml}
\alias{mvmeta.ml}
\alias{mvmeta.reml}
\alias{mvmeta.fixed}


\title{ Fitting Functions for mvmeta Models }

\description{
These functions defines alternative fitting methods for multivariate and univariate meta-analysis and meta-regression. They are meant to be used internally and not directly run by the users.
}

\usage{
mvmeta.fixed(Xlist, ylist, Slist, nalist, k, m, nall, control)

mvmeta.ml(Xlist, ylist, Slist, nalist, k, m, nall, control)

mvmeta.reml(Xlist, ylist, Slist, nalist, k, m, nall, control)
}

\arguments{
Assuming a meta-analysis or meta-regression based on \eqn{m} studies, \eqn{k} outcomes and \eqn{p} predictors:

  \item{Xlist }{ a \eqn{m}-dimensional list of study-specific design matrices for the fixed-effects part of the model. Rows corresponding to missing outcomes have been excluded.}
  \item{ylist }{ a \eqn{m}-dimensional list of study-specific vectors of estimated outcomes. Entries corresponding to missing outcomes have been excluded.}
  \item{Slist }{ a \eqn{m}-dimensional list of within-study (co)variance matrices of estimated outcomes. Rows and columns corresponding to missing outcomes have been excluded.}
  \item{nalist }{ a \eqn{m}-dimensional list of \eqn{k}-dimensional study-specific logical vectors, identifying missing outcomes.}
  \item{k, m, nall }{ numeric scalars: number of outcomes, number of studies included in estimation (equal to the length of lists above), number of observations (excluding missing).}
  \item{control }{ list of parameters for controlling the fitting process, usually provided by \code{\link{glm.control}}.}
}

\details{
Fixed-effects meta-analysis and meta-regression models are fitted by \code{mvmeta.fitted}, using a generalized least square algorithm implemented in the internal function \code{.gls}. The estimation procedure focuses on the \eqn{kp} fixed-effects coefficients, while assuming \code{Psi} as a 0 matrix.

Random-effects models are fitted through profiled (restricted) maximum likelihood by the functions \code{mvmeta.ml} and \code{mvmeta.reml}. They perform a hybrid estimation procedure based first on few runs of iterative generalized least square algorithm, implemented in \code{\link{mvmeta.igls}}, and then quasi-Newton iterations until convergence, provided by the general-purpose optimizer \code{\link{optim}} based on specific \code{\link[=mvmeta.ml.fn]{likelihood functions}}. The estimation focuses on the \eqn{kp} fixed-effects coefficients and the \eqn{k(k+1)/2} variance components for the random part of the model, corresponding to the lower triangular entries of the between-study (co)variance matrix \code{Psi}.

As explained above, the arguments of these functions are prepared by \code{\link{mvmeta.fit}} from its arguments \code{X}, \code{y} and \code{S}, plus a control list passed by \code{\link{glm.control}}. The list structure, although requiring more elaborate coding, is computationally more efficient, as it avoids the specification of sparse block-diagonal matrices, especially for meta-analysis involving a large number of studies.

In the absence of missing values, each \eqn{k \times kp}{k x kp} matrix in \code{Xlist} is produced by a Kronecker product between an identity matrix of dimension \eqn{k} and the transpose of the \eqn{p}-dimensional vector of study-level predictors, corresponding to the related row of \code{X}. The \eqn{k}-dimensional vectors in \code{ylist} are simply the entries of the \eqn{m}-dimensional vector (for univariate models) or rows of the \eqn{m \times k}{m x k} matrix (for multivariate models) of outcomes \code{y}. The within-study (co)variance matrices in \code{Slit} are produced by expanding the rows of the \eqn{k \times k(k+1)/2}{k x k(k+1)/2} matrix \code{S}, stored in vectorized form.  In the presence of missing outcomes in \code{y}, the related parts of the each component of these lists are excluded, as no missing values are accepted in these fitting functions and related \code{\link[=mvmeta.ml.fn]{likelihood functions}}.
}

\value{
These functions return an intermediate list object, whose components are then processed by \code{\link{mvmeta.fit}}. Other components are added later through \code{\link{mvmeta}} to finalize an object of class \code{"mvmeta"} .
}

\references{
Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. Epub ahead of print (doi 10.1002/sim.5471). [Freely available \href{http://onlinelibrary.wiley.com/doi/10.1002/sim.5471/pdf}{\bold{here}}].

Goldstein H (1986). Multilevel mixed linear model analysis using iterative generalized least squares. \emph{Biometrika}. \bold{73}(1):43.

Pinheiro JC and Bates DM (2000). Mixed-Effects Models in S and S-PLUS. New York, Springer Verlag.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
As stated earlier, these functions are called internally by \code{\link{mvmeta.fit}}, and are not meant to be used directly. In particular, their code does not contain any check on the arguments provided, which are expected in specific formats. They are however exported in the namespace and documented for completeness.
}

\seealso{
See \code{\link{mvmeta.fit}} and the \code{\link[=mvmeta.ml.fn]{likelihood functions}} for additional info on the fitting procedures. See the function \code{\link{getAnywhere}} to access the code of the internal functions, such as \code{.gls}. See \code{\link{mvmeta-package}} for an overview of the package and modelling framework.
}


\keyword{models}
\keyword{regression}
\keyword{multivariate}
