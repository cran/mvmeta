\name{mvmeta.ml}
\alias{mvmeta.ml}
\alias{mvmeta.reml}
\alias{mvmeta.fixed}
\alias{mvmeta.mm}


\title{ Fitting Functions for mvmeta Models }

\description{
These functions defines alternative fitting methods for multivariate and univariate meta-analysis and meta-regression. They are meant to be used internally and not directly run by the users.
}

\usage{
mvmeta.fixed(Xlist, ylist, Slist, nalist, k, m, p, nall, control)

mvmeta.ml(Xlist, ylist, Slist, nalist, k, m, p, nall, control)

mvmeta.reml(Xlist, ylist, Slist, nalist, k, m, p, nall, control)

mvmeta.mm(Xlist, ylist, Slist, nalist, k, m, p, nall, control)
}

\arguments{
Assuming a meta-analysis or meta-regression based on \eqn{m} studies, \eqn{k} outcomes and \eqn{p} predictors:

  \item{Xlist }{ a \eqn{m}-dimensional list of study-specific design matrices for the fixed-effects part of the model. Rows corresponding to missing outcomes have been excluded.}
  \item{ylist }{ a \eqn{m}-dimensional list of study-specific vectors of estimated outcomes. Entries corresponding to missing outcomes have been excluded.}
  \item{Slist }{ a \eqn{m}-dimensional list of within-study (co)variance matrices of estimated outcomes. Rows and columns corresponding to missing outcomes have been excluded.}
  \item{nalist }{ a \eqn{m}-dimensional list of \eqn{k}-dimensional study-specific logical vectors, identifying missing outcomes.}
  \item{k, m, p, nall }{ numeric scalars: number of outcomes, number of studies included in estimation (equal to the length of lists above), number of predictors (including the intercept), number of observations (excluding missing).}
  \item{control }{ list of parameters for controlling the fitting process, usually internally set to default values by \code{\link{mvmeta.control}}.}
}

\details{
Meta-analytical models divide into fixed and random-effects models, the latter fitted with different procedures based here on likelihood or method of moments. In fixed-effects models, the estimation only involves \eqn{kp} coefficients, while in random-effects models also \eqn{k(k+1)/2} variance components are estimated, corresponding to the lower triangular entries of the between-study (co)variance matrix \code{Psi}. 

Each fitting function performs a different estimation method. Specifically:

\itemize{
  \item \code{mvmeta.fitted} fits fixed-effects meta-analysis and meta-regression. The between-study (co)variance matrix is set to zero, so the marginal (co)variance matrix, composed only by elements of the within-study component, is assumed as completely known. The routine is based on a standard generalized least square (GLS) algorithm implemented in the internal function \code{.gls}. Similarly to the likelihood method below, the computation involves Cholesky and and QR decompositions for computational stability and efficiency (see the help page for the \code{\link[=mvmeta.ml.fn]{likelihood functions}}). The method is described in details in Gasparrini and collaborators (2012) (see references below).
  \item \code{mvmeta.ml} and \code{mvmeta.reml} fit random-effects meta-analysis and meta-regression through maximum likelihood (ML) and restricted maximum likelihood (REML). They perform a hybrid estimation procedure based first on few runs of iterative generalized least square algorithm and then quasi-Newton iterations, using specific \code{\link[=mvmeta.ml.fn]{likelihood functions}}, until convergence. The estimation algorithm adopts a profiled (or concentrated) approach, that is expressed only in terms of the variance components. Cholesky and and QR decompositions are used for computational stability and efficiency, and for assuring the positive-definiteness of the estimated between-study (co)variance matrix. See the help page for the \code{\link[=mvmeta.ml.fn]{likelihood functions}} for further details. The method is described in details in Gasparrini and collaborators (2012) (see references below).
  \item \code{mvmeta.mm} fits random-effects meta-analysis and meta-regression through the method of moments. This approach represents the multivariate extension of the traditional estimator used in meta-analysis, and simplifies to the standard method in the univariate case. This non-iterative routine assures the positive semi-definiteness of the estimated between-study (co)variance matrix by setting its negative eigenvalues to zero. The specific method of moment estimator used here is described in Chen and collaborators (2012) and Jackson and collaborators (2013) (see references below).
}

Some parameters of the fitting procedures may be determined through \code{\link{mvmeta.control}}. 

As explained above, the arguments of these functions are prepared by \code{\link{mvmeta.fit}} from its arguments \code{X}, \code{y} and \code{S}, plus a control list passed by \code{\link{glm.control}}. The list structure, although requiring more elaborate coding, is computationally more efficient, as it avoids the specification of sparse block-diagonal matrices, especially for meta-analysis involving a large number of studies.

In the absence of missing values, each \eqn{k \times kp}{k x kp} matrix in \code{Xlist} is produced by a Kronecker product between an identity matrix of dimension \eqn{k} and the transpose of the \eqn{p}-dimensional vector of study-level predictors. The \eqn{k}-dimensional vectors in \code{ylist} are simply the entries of the \eqn{m}-dimensional vector (for univariate models) or rows of the \eqn{m \times k}{m x k} matrix (for multivariate models) of outcomes \code{y}. The within-study (co)variance matrices in \code{Slist} are produced by expanding the rows of the \eqn{k \times k(k+1)/2}{k x k(k+1)/2} matrix \code{S}, stored in vectorized form.  In the presence of missing outcomes in \code{y}, the related parts of the each component of these lists are excluded, as no missing values are accepted in these fitting functions and related \code{\link[=mvmeta.ml.fn]{likelihood functions}}.
}

\value{
These functions return an intermediate list object, whose components are then processed by \code{\link{mvmeta.fit}}. Other components are added later through \code{\link{mvmeta}} to finalize an object of class \code{"mvmeta"} .
}

\references{
Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/statmed2012.html}{\bold{here}}].

Jackson D, White IR, Riley RD (2013). A matrix based method of moments for fitting the multivariate random effects model for meta-analysis and meta-regression. \emph{Biometrical Journal}. To be published.

Chen H, Manning AK, Dupuis J (2012). A method of moments estimator for random effect multivariate meta-analysis. \emph{Biometrics}. \bold{68}(4):1278-1284.

Goldstein H (1986). Multilevel mixed linear model analysis using iterative generalized least squares. \emph{Biometrika}. \bold{73}(1):43.

Pinheiro JC and Bates DM (2000). Mixed-Effects Models in S and S-PLUS. New York, Springer Verlag.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
As stated earlier, these functions are called internally by \code{\link{mvmeta.fit}}, and are not meant to be used directly. In particular, their code does not contain any check on the arguments provided, which are expected in specific formats. They are however exported in the namespace and documented for completeness.
}

\seealso{
See \code{\link{mvmeta.fit}} and the \code{\link[=mvmeta.ml.fn]{likelihood functions}} for additional info on the fitting procedures. See \code{\link{mvmeta.control}} to determine specific parameters of the fitting procedures. See the function \code{\link{getAnywhere}} to access the code of the internal functions, such as \code{.gls}. See \code{\link{mvmeta-package}} for an overview of the package and modelling framework.
}


\keyword{models}
\keyword{regression}
\keyword{multivariate}
