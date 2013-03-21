\name{mvmeta.control}
\alias{mvmeta.control}


\title{ Ancillary Parameters for Controlling the Fit in mvmeta Models }

\description{
This internal function sets the parameter options used for fitting meta-analytical models, commonly to pre-specified default values. It is usually internally called by \code{\link{mvmeta.fit}}.
}

\usage{
mvmeta.control(optim=list(), showiter=FALSE, maxiter=100, igls.iter=10,
  reltol=sqrt(.Machine$double.eps),vc.adj=TRUE,
  set.negeigen=sqrt(.Machine$double.eps))
}

\arguments{

  \item{optim }{ List of parameters passed to the \code{control} argument of the function \code{\link{optim}}, which performs the quasi-Newton optimization in likelihood-based random-effects models. See \code{\link{optim}}.}
  \item{showiter }{ Logical. If \code{TRUE}, the progress of iterative optimization is shown.}
  \item{maxiter }{ Positive interger value. Maximum number of iterations in methods involving optimization procedures.}
  \item{igls.iter }{ Number of iteration of the iterative generalized least square algorithm to be run in the hybrid optimization procedure used in likelihood-based random-effects model.}
  \item{reltol }{ Relative convergence tolerance in methods involving optimization procedures. The algorithm stops if it is unable to reduce the value by a factor of \code{reltol * (abs(val) + reltol)} at a step.}
  \item{vc.adj }{ Logical. If \code{TRUE}, an adjustement to the way the marginal variance part is computed in the variance components estimator is applied. See \code{\link{mvmeta.vc}}.}
  \item{set.negeigen }{ Positive value. Value to which negative eigenvalues are to be set in estimators where such method is used to force positive semi-definiteness of the estimated between-study (co)variance matrix.}
}

\details{
The control argument of \code{\link{mvmeta}} is by default passed to \code{\link{mvmeta.fit}}, which uses its elements as arguments of \code{mvmeta.control}.

Many arguments refer to specific fitting procedures. Refer to the help page of the related estimator for details.

The function is expected to be extended and/or modified at every releases of the package \pkg{mvmeta}.
}

\value{
A list with components named as the arguments.
}

\references{
Gasparrini A, Armstrong B, Kenward MG (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. \bold{31}(29):3821--3839. [Freely available \href{http://www.ag-myresearch.com/statmed2012.html}{\bold{here}}].
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\seealso{ 
See \code{\link{mvmeta}}. See also \code{\link{glm.control}}. See the help pages of the related fitting functions for details on each parameter. See \code{\link{mvmeta-package}} for an overview of this modelling framework.
}

\examples{
# RUN THE MODEL
model <- mvmeta(cbind(PD,AL)~pubyear,S=berkey98[5:7],data=berkey98)

# PRINT THE ITERATIONS (SEE ?optim)
model <- mvmeta(cbind(PD,AL)~pubyear,S=berkey98[5:7],data=berkey98,
  control=list(showiter=TRUE,igls.iter=20))
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
