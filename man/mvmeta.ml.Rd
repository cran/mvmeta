\name{mvmeta.ml}
\alias{mvmeta.ml}
\alias{mvmeta.reml}
\alias{mvmeta.ml.grad}
\alias{mvmeta.reml.grad}
\alias{mvmeta.igls}


\title{ Internal Functions for Package mvmeta }

\description{
These functions are used in algorithms for fitting random-effects multivariate meta-analysis and meta-regression through (restricted) maximum likelihood. They are called internally by \code{\link{mvmeta}}, and are not meant to be directly run by the users (see Note below).
}

\usage{

mvmeta.ml(par, ylist, Slist, kXlist, nalist, nobs, k)
mvmeta.reml(par, ylist, Slist, kXlist, nalist, nobs, k)

mvmeta.ml.grad(par, ylist, Slist, kXlist, nalist, nobs, k)
mvmeta.reml.grad(par, ylist, Slist, kXlist, nalist, nobs, k)

mvmeta.igls(Psi, ylist, Slist, kXlist, nalist, k, m, p)

}

\arguments{

  \item{par }{ numerical vector representing the upper triangle entries of the Cholesky decomposition of \code{Psi}, taken by row.}
  \item{Psi }{ matrix representing the current estimate of the between-study (co)variance matrix of the outcome parameters.}
  \item{ylist }{ a list of study-specific vectors of estimated outcome parameters. Missing values have been excluded.}
  \item{Slist }{ a list of study-specific (co)variance matrices of estimated outcome parameters. Missing values have been excluded.}
  \item{kXlist }{ a list of study-specific design matrices for the fixed part of the model, as returned by \code{\link{kXlistmk}}.}
  \item{nalist }{ a list of \eqn{k}-length study-specific logical vectors, \code{TRUE} when the outcome parameter is not missing.}
  \item{nobs, k, m, p }{ numeric scalars: number of observations (excluding missing), number of outcome parameters, number of studies included in estimation (equal to the length of lists above), number of coefficients for each outcome parameter.}

}

\details{
Random-effects models are fitted through (restricted) maximum likelihood using a quasi-Newton iterative algorithm, in order to estimate the coefficients \eqn{beta} and the variance components of \eqn{Psi}. The procedure is based on a profiled (or concentrated) likelihood, expressing the (restricted) likelihood only in terms of the components of \eqn{Psi}, while the conditional estimate of \eqn{beta}, obtained through generalized least square, is plugged-in at every iteration. This approach reduces the number of parameters in the optimization routine, especially in meta-regression models.

The optimization routine is performed by \code{\link{optim}}, which calls \code{mvmeta.ml}-\code{mvmeta.ml.grad} or \code{mvmeta.reml}-\code{mvmeta.reml.grad}. These functions provide the (restricted) likelihood functions and their functions of partial derivatives with respect to the model parameters. Given that \code{\link{optim}} minimizes the objective, these values are returned with the minus sign.
Although efficient, convergence of Newton methods is highly dependent on optimal starting values of the parameters, which are computed by few iterations of an iterative generalized least square algorithm provided by the function \code{mvmeta.igls}. Initial value of \eqn{Psi} for \code{mvmeta.igls} is given by a diagonal matrix with 0.001 entries.

The parameters of \eqn{Psi} are expressed in terms of its Cholesky decomposition, in order to assure the positive-definiteness of the estimated between-study (co)variance matrix. (Restricted) likelihood and derivatives are computed as sums of study-level contributions, in order to avoid the specification of high-dimensional sparse matrices and reduce the computation burden. A second Cholesky decomposition of  \eqn{Sigma = S +Psi} is performed in order to speed up the computation of \eqn{beta} through generalized least-square at each iteration, and the calculation of matrix inverses and determinants. More details are provided in the references below.
}

\value{
\code{mvmeta.ml} and \code{mvmeta.reml} return the value of minus the (restricted) log-likelihood for a given set of parameters in \code{par}. \code{mvmeta.ml.grad} and \code{mvmeta.reml.grad} return instead the related vector of partial derivatives. \code{mvmeta.igls} returns an updated estimate of \eqn{Psi} given its initial value or the value at previous iteration.
}

\references{
Goldstein, H. (1986). "Multilevel mixed linear model analysis using iterative generalized least squares." \emph{Biometrika}. \bold{73}(1):43.

Lindstrom, M. J. and D. M. Bates (1988). "Newton-Raphson and EM algorithms for linear mixed-effects models for repeated-measures data." \emph{Journal of the American Statistical Association}. \bold{83}(404):1014-1022.

Pinheiro, J. C. and D. M. Bates (2000). Mixed-Effects Models in S and S-PLUS. New York, Springer Verlag.
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{As stated earlier, these functions are called internally within the \pkg{mvmeta} package, and are not meant to be used directly. In particular, their code does not contain any check on the arguments provided, which are expected in specific formats. Although not exported (included in the namespace), they are documented in order to facilitate inspection and intepretation of the package code. Their code is reported below.
}

\seealso{
See \code{\link{mvmeta-package}} and \code{\link{mvmeta}} for an overview of the package and modelling framework.
}

\examples{
mvmeta.ml <- function(par, ylist, Slist, kXlist, nalist, nobs, k) {
	# PARAMETERIZATION OF theta AS THE LOWER TRIANGULAR COMPONENTS OF
	#	THE SQUARE ROOT OF Psi, THROUGH CHOLESKY-DECOMPOSITION
	Psi <- matrix(0,k,k)
	Psi[lower.tri(Psi,diag=TRUE)] <- par
	Psi <- tcrossprod(Psi)

	# COMPUTE beta BY GLS
	Sigmalist <- mapply(function(S,na) S+Psi[na,na,drop=FALSE],
		Slist,nalist,SIMPLIFY=FALSE)
	Ulist <- lapply(Sigmalist,chol)
	invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
	invtUXlist <- mapply(function(invU,kX) crossprod(invU,kX),
		invUlist,kXlist,SIMPLIFY=FALSE)
	invtUylist <- mapply(function(invU,y) crossprod(invU,y),	
		invUlist,ylist,SIMPLIFY=FALSE)
	invtUX <- rbindlist(invtUXlist)
	invtUy <- rbindlist(invtUylist)
	beta <- as.numeric(qr.solve(invtUX,invtUy))
	# LIKELIHOOD FUNCTION
	# CONSTANT PART
	const <- -0.5*nobs*log(2*pi)
	# I GUESS IN STATA:
	#const <- -0.5*length(ylist)*ncol(Psi)*log(2*pi)
	# RESIDUAL COMPONENT
	res <- -0.5*crossprod(invtUy-invtUX\%*\%beta)
	# DETERMINANT COMPONENT
	det <- -sum(sapply(Ulist,function(U) sum(log(diag(U)))))
	logLik <- const + det + res
	# RETURN MINUS THE VALUE: optim MINIMIZES
	return(-logLik)
}

mvmeta.reml <- function(par, ylist, Slist, kXlist, nalist, nobs, k) {
	# PARAMETERIZATION OF theta AS THE LOWER TRIANGULAR COMPONENTS OF
	#	THE SQUARE ROOT OF Psi, THROUGH CHOLESKY-DECOMPOSITION
	Psi <- matrix(0,k,k)
	Psi[lower.tri(Psi,diag=TRUE)] <- par
	Psi <- tcrossprod(Psi)
	# GET THE ESTIMATE OF beta CONDITIONAL ON Psi_theta
	Sigmalist <- mapply(function(S,na) S+Psi[na,na,drop=FALSE],
		Slist,nalist,SIMPLIFY=FALSE)
	Ulist <- lapply(Sigmalist,chol)
	invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
	invtUXlist <- mapply(function(invU,kX) crossprod(invU,kX),
		invUlist,kXlist,SIMPLIFY=FALSE)
	invtUylist <- mapply(function(invU,y) crossprod(invU,y),
		invUlist,ylist,SIMPLIFY=FALSE)
	beta <- qr.solve(rbindlist(invtUXlist),rbindlist(invtUylist))
	# LIKELIHOOD FUNCTION
	# CONSTANT PART
	const <- -0.5*(nobs-length(beta))*log(2*pi)
	# I GUESS IN STATA:
	#const <- -0.5*(length(ylist)-1)*ncol(Psi)*log(2*pi)
	# RESIDUAL COMPONENT
	res <- -0.5*sum(mapply(function(invtUy,invtUX) {
		crossprod(invtUy-invtUX\%*\%beta)},invtUylist,invtUXlist))
	# DETERMINANT COMPONENTS
	det1 <- -sum(sapply(Ulist,function(U) sum(log(diag(U)))))
	tXMXtot <- sumlist(lapply(invtUXlist,function(x)crossprod(x)))
	det2 <- -sum(log(diag(chol(tXMXtot))))
	logLik <- const + det1 + det2 + res
	# RETURN MINUS THE VALUE: optim MINIMIZES
	return(-logLik)
}

mvmeta.ml.grad <- function(par, ylist, Slist, kXlist, nalist, nobs, k) {
	# RETRIEVE THE UPPER TRIANGULAR CHOLESKY MATRIX AND Psi
	L <- matrix(0,k,k)
	L[lower.tri(L,diag=TRUE)] <- par
	U <- t(L)
	Psi <- crossprod(U)
	# GET THE ESTIMATE OF beta CONDITIONAL ON Psi_theta
	Sigmalist <- mapply(function(S,na) S+Psi[na,na,drop=FALSE],
		Slist,nalist,SIMPLIFY=FALSE)
	Ulist <- lapply(Sigmalist,chol)
	invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
	invtUXlist <- mapply(function(invU,kX) crossprod(invU,kX),
		invUlist,kXlist,SIMPLIFY=FALSE)
	invtUylist <- mapply(function(invU,y) crossprod(invU,y),	
		invUlist,ylist,SIMPLIFY=FALSE)
	invtUX <- rbindlist(invtUXlist)
	invtUy <- rbindlist(invtUylist)
	beta <- as.numeric(qr.solve(invtUX,invtUy))
	# COMPUTE THE MATRIX DERIVATIVES OF EACH PARAMETER
	Dlist <- lapply(seq(length(par)),gradchol,k,U)
	# FUNCTION TO COMPUTE THE GRADIENT COMPONENTS FOR A GIVEN PARAMETER
	fgrad <- function(D,y,kX,invU,na) {
		sum(mapply(function(y,kX,invU,na) {
			invSigma <- tcrossprod(invU)
			res <- as.numeric(y-kX\%*\%beta)
			A <- crossprod(res,invSigma)\%*\%D[na,na]%*%invSigma\%*\%res
			B <- sum(diag(invSigma\%*\%D[na,na]))
			return(as.numeric(0.5*(A-B)))},y,kX,invU,na))
	}
	# REPEAT FOR ALL THE PARAMETERS 
	grad <- sapply(Dlist,fgrad,y=ylist,kX=kXlist,invU=invUlist,nalist)
	return(-grad)
}

mvmeta.reml.grad <- function(par, ylist, Slist, kXlist, nalist, nobs, k) {
	# RETRIEVE THE UPPER TRIANGULAR CHOLESKY MATRIX AND Psi
	L <- matrix(0,k,k)
	L[lower.tri(L,diag=TRUE)] <- par
	U <- t(L)
	Psi <- crossprod(U)
	# GET THE ESTIMATE OF beta CONDITIONAL ON Psi_theta
	Sigmalist <- mapply(function(S,na) S+Psi[na,na,drop=FALSE],
		Slist,nalist,SIMPLIFY=FALSE)
	Ulist <- lapply(Sigmalist,chol)
	invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
	invtUXlist <- mapply(function(invU,kX) crossprod(invU,kX),
		invUlist,kXlist,SIMPLIFY=FALSE)
	invtUylist <- mapply(function(invU,y) crossprod(invU,y),
		invUlist,ylist,SIMPLIFY=FALSE)
	beta <- qr.solve(rbindlist(invtUXlist),rbindlist(invtUylist))

	tXMXtot <- sumlist(lapply(invtUXlist,function(x)crossprod(x)))
	invtXMXtot <- chol2inv(chol(tXMXtot))
	# COMPUTE THE MATRIX DERIVATIVES OF EACH PARAMETER
	Dlist <- lapply(seq(length(par)),gradchol,k,U)
	# FUNCTION TO COMPUTE THE GRADIENT COMPONENTS FOR A GIVEN PARAMETERS
	fgrad <- function(D,y,kX,invU,na) {
		sum(mapply(function(y,kX,invU,na) {
			invSigma <- tcrossprod(invU)
			MDM <- invSigma\%*\%D[na,na]\%*\%invSigma
			res <- as.numeric(y-kX\%*\%beta)
			A <- crossprod(res,MDM)\%*\%res
			B <- sum(diag(invSigma\%*\%D[na,na]))
			C <- sum(diag(invtXMXtot\%*\%crossprod(kX,MDM)\%*\%kX))
			return(as.numeric(0.5*(A-B+C)))},y,kX,invU,na))
	}
	# REPEAT FOR ALL THE PARAMETERS 
	grad <- sapply(Dlist,fgrad,y=ylist,kX=kXlist,invU=invUlist,nalist)
	return(-grad)
}

mvmeta.igls <- function(Psi, ylist, Slist, kXlist, nalist, k, m, p) {
	# COMPUTE beta BY GLS
	Sigmalist <- mapply(function(S,na) S+Psi[na,na,drop=FALSE],
		Slist,nalist,SIMPLIFY=FALSE)
	Ulist <- lapply(Sigmalist,chol)
	invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
	invtUXlist <- mapply(function(invU,kX) crossprod(invU,kX),
		invUlist,kXlist,SIMPLIFY=FALSE)
	invtUylist <- mapply(function(invU,y) crossprod(invU,y),	
		invUlist,ylist,SIMPLIFY=FALSE)
	invtUX <- rbindlist(invtUXlist)
	invtUy <- rbindlist(invtUylist)
	beta <- as.numeric(qr.solve(invtUX,invtUy))
	# CREATE A MATRIX WITH INDICATOR OF (CO)VAR COMPONENTS
	npar <- k*(k+1)/2
	indMat <- xpndMat(seq(npar))
	# CREATE THE TRANSFORMED OBJECTS IN IGLS STRUCTURE
	# EXPANDED (CO)VARIANCE MATRIX
	Sigmalist <- mapply(function(na,S) {
		V <- Psi[na,na,drop=FALSE]+S
		return(V\%x\%V)},nalist,Slist,SIMPLIFY=FALSE)
	# RESPONSE VECTORS WITH RESIDUALS MINUS THE WITHIN (CO)VARIANCE
	#	COMPONENTS, CONSIDERED FIXED
	flist <- mapply(function(y,S,kX) {
		return(as.numeric(tcrossprod(y-kX\%*\%beta))-as.numeric(S))},
		ylist,Slist,kXlist,SIMPLIFY=FALSE)
	# DESIGN MATRIX MAPPING THE PARAMETERS TO BE ESTIMATED
	#	IT AUTOMATICALLY CREATES 0 COLUMNS FOR MISSING OBSERVATIONS
	Zlist <- lapply(nalist,function(na) {
		z <- as.numeric(indMat[na,na,drop=FALSE])
		Z <- lapply(seq(npar),function(x) as.numeric(z==x))
		return(cbindlist(Z))})
	# CREATE TRANFORMED OBJECTS FOR WEIGHTED LEAST-SQUARE THROUGH CHOLESKY
	Ulist <- lapply(Sigmalist,chol)
	invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
	invtUZlist <- mapply(function(invU,Z) crossprod(invU,Z),
		invUlist,Zlist,SIMPLIFY=FALSE)
	invtUflist <- mapply(function(invU,f) crossprod(invU,f),
		invUlist,flist,SIMPLIFY=FALSE)
	invtUZ <- rbindlist(invtUZlist)
	invtUf <- rbindlist(invtUflist)
	# ESTIMATE THE COMPONENTS
	theta <- as.numeric(qr.solve(invtUZ,invtUf))
	Psi <- xpndMat(theta)
	# FORCING POSITIVE-DEFINITENESS
	Psi <- expMat(logMat(Psi,positive=TRUE))
	return(Psi)
}

}

\keyword{internal}
\keyword{models}
\keyword{regression}
\keyword{multivariate}
