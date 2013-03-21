###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`.gls` <-
function(Xlist, ylist, Slist, nalist, Psi, onlycoef=TRUE)  {
#
################################################################################
# FUNCTION TO COMPUTE THE GLS ESTIMATE + OPTIONAL INTERMEDIATE PRODUCTS
#
  Sigmalist <- mapply(function(S,na) S+Psi[!na,!na,drop=FALSE],
    Slist,nalist,SIMPLIFY=FALSE)
  Ulist <- lapply(Sigmalist,chol)
  invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
  invtUXlist <- mapply(function(invU,X) crossprod(invU,X),
    invUlist,Xlist,SIMPLIFY=FALSE)
  invtUylist <- mapply(function(invU,y) crossprod(invU,y),
    invUlist,ylist,SIMPLIFY=FALSE)
  invtUX <- do.call("rbind",invtUXlist)
  invtUy <- do.call("rbind",invtUylist)
  coef <- as.numeric(qr.solve(invtUX,invtUy))
#
  if(onlycoef) return(coef)
  return(list(coef=coef,Sigmalist=Sigmalist,Ulist=Ulist,invUlist=invUlist,
    invtUXlist=invtUXlist,invtUX=invtUX,invtUy=invtUy))
}
#
###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`.gradchol.ml` <-
function(par,U,ind1,ind2,invSigmalist,reslist,nalist,k,m) {
#
################################################################################
# FUNCTION TO COMPUTE THE MATRIX DERIVATIVES IN TERMS OF PARAMETERS OF
# THE CHOLESKY DECOMPOSITION
#
  grad <- sapply(seq(length(par)),function(i) {
    # COMPUTE THE DERIVATIVE OF Psi IN TERMS OF ITS CHOLESKY DECOMPOSITION U
    A <- B <- C <- diag(0,k)
    A[ind2[i],] <- B[,ind2[i]] <- U[ind1[i],]
    C[ind2[i],] <- C[,ind2[i]] <- 1
    D <- C*A+C*B
    # COMPUTE THE GRADIENT
    gr <- sum(mapply(function(invSigma,res,na) {
      E <- crossprod(res,invSigma)%*%D[!na,!na]%*%invSigma%*%res
      F <- sum(diag(invSigma%*%D[!na,!na]))
      return(as.numeric(0.5*(E-F)))},invSigmalist,reslist,nalist))
  })
#
  return(grad)
}
###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`.gradchol.reml` <-
function(par,U,invtXWXtot,ind1,ind2,Xlist,invSigmalist,reslist,nalist,k,m) {
#
################################################################################
# FUNCTION TO COMPUTE THE MATRIX DERIVATIVES IN TERMS OF PARAMETERS OF
# THE CHOLESKY DECOMPOSITION
#
  grad <- sapply(seq(length(par)),function(i) {
    # COMPUTE THE DERIVATIVE OF Psi IN TERMS OF ITS CHOLESKY DECOMPOSITION U
    A <- B <- C <- diag(0,k)
    A[ind2[i],] <- B[,ind2[i]] <- U[ind1[i],]
    C[ind2[i],] <- C[,ind2[i]] <- 1
    D <- C*A+C*B
    # COMPUTE THE GRADIENT
    gr <- sum(mapply(function(X,invSigma,res,na) {
      E <- invSigma%*%D[!na,!na]%*%invSigma
      F <- crossprod(res,E)%*%res
      G <- sum(diag(invSigma%*%D[!na,!na]))
      H <- sum(diag(invtXWXtot%*%crossprod(X,E)%*%X))
      return(as.numeric(0.5*(F-G+H)))},Xlist,invSigmalist,reslist,nalist))
  })
#
  return(grad)
}
###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`.mkS` <- 
function(S, y, narm=NULL, subset=NULL) {
#
################################################################################
# TRANSFORM S IN A MATRIX OF VECTORIZED VCOV MATRICES
#
  # DIMENSIONS
  k <- dim(y)[2]
  m <- dim(y)[1]
#
  # MESSAGE
  mes <- "incorrect dimensions for 'S'"
#
  # IF DATAFRAME
  if(is.data.frame(S)) S <- as.matrix(S)
#
  # IF NUMERIC
  if(is.numeric(S)) {
    # IF JUST A VECTOR, TRANSFORM IN A MATRIX
    if(is.null(dim(S))) S <- as.matrix(S)
    # IF AN ARRAY, TRANSFORM IN A MATRIX
    if(length(dim(S))==3L) S <- t(apply(S,3,vechMat))
    # FINALLY, IF A MATRIX
    if(!is.null(subset)) S <- S[subset,,drop=FALSE]
    if(!is.null(narm)) S <- S[-narm,,drop=FALSE]
    if(dim(S)[1]!=m || dim(S)[2]!=k*(k+1)/2) stop(mes)    
  }
#
  # IF A LIST
  if(is.list(S)) {
    if(!is.null(subset)) S <- S[subset]
    if(!is.null(narm)) S <- S[-narm]
    if(length(S)!=m) stop(mes)
    if(any(sapply(S,dim)!=k)) stop(mes)
    S <- if(k==1L) as.matrix(sapply(S,vechMat)) else t(sapply(S,vechMat))
  }
#
  # NAMES
  rownames(S) <- rownames(y)
  nk <- colnames(y)
  colnames(S) <- if(!is.null(nk)) vechMat(outer(nk,nk,paste,sep=".")) else NULL
#
  return(S)
}
###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`.mvsim` <- 
function(nsim=1, mu, Sigma, posdeftol=sqrt(.Machine$double.eps), drop=TRUE) {
#
################################################################################
#
  # DIMENSIONS
  k <- length(mu)
  eigen <- eigen(Sigma,symmetric=TRUE)
  # EIGENVALUE DECOMPOSITION
  ev <- eigen$values
#
  # CHECK POSITIVE DEFINITENESS
  if(any(ev < -posdeftol*abs(max(ev)))) stop("lack of positive definiteness")
#
  # SIMULATE
  X <- matrix(rnorm(k*nsim),nsim)
  sim <- drop(mu) + eigen$vectors %*% diag(sqrt(pmax(ev,0)),k) %*% t(X)
#
  if(drop) return(drop(t(sim))) else return(t(sim)) 
}
###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`.onAttach` <-
function(lib, pkg) {
#
################################################################################
#
  meta <- packageDescription("mvmeta")
  attachmsg <- paste("This is mvmeta ",meta$Version,
    ". For an overview type: help('mvmeta-package').",sep="")
  packageStartupMessage(attachmsg, domain = NULL, appendLF = TRUE)
}
###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`.sumlist` <-
function(list) {
#
################################################################################
#
  res <- 0
  for(i in seq(list)) res <- res + list[[i]]
  return(res)
}
###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`.fbtr` <- 
function(A,k) {
#
################################################################################
#
  btrA <- 0
  for(i in seq(dim(A)[1]/k)) {
    ind <- (i-1)*k+1:k
    btrA <- btrA + A[ind,ind]
  }
  return(btrA)
}
