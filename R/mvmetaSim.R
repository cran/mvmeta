###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2013
#
`mvmetaSim` <- 
function(y, S, Psi, sd, cor, nsim=1, seed=NULL, posdeftol) {
#
################################################################################
#
  if(missing(posdeftol)) posdeftol <- sqrt(.Machine$double.eps)
#
  # DEFINE THE SEED (FROM simulate.lm)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
#
  # PREPARE AND CHECK Psi, THEN SET THE DIMENSION
  if(missing(Psi)) {
    if(missing(sd)||missing(cor)) stop("'Psi' or 'sd'-'cor' must be provided")
    R <- diag(1,length(sd))
    R[lower.tri(R)] <- if(is.matrix(cor)) cor[lower.tri(cor)] else cor
    R[upper.tri(R)] <- t(R)[upper.tri(t(R))]
    if(any(R^2>1)) stop("correlations must be between -1 and 1")
    Psi <- diag(sd)%*%R%*%diag(sd)
  }
  if(!is.matrix(Psi)) Psi <- xpndMat(Psi)
  dim <- dim(Psi)
  if(dim[1]!=dim[2]) stop("'Psi' is not a square matrix")
  if(any(is.na(Psi))) stop("missing values not allowed in 'Psi'")
  k <- dim[1]
#
  # PREPARE AND CHECK y
  if(!is.matrix(y)) y <- as.matrix(y)
  if(ncol(y)!=k) stop("Dimensions of 'y' and 'Psi' not consistent")
  if(any(is.na(y))) stop("missing values not allowed in 'y'")
#
  # PREPARE AND CHECK S
  S <- .mkS(S,y)
  if(any(is.na(S))) stop("missing values not allowed in 'S'")
#
  # SAMPLE THE RESPONSES
  # FOR EFFICIENCY, IT SAMPLES SEVERAL OUTCOMES FROM THE SAME MEAN AND
  #   THEN RE-ARRANGE THEM
  sim <- do.call("cbind",lapply(seq(nrow(y)), function(i) {
    .mvsim(nsim,y[i,],Sigma=xpndMat(S[i,])+Psi,posdeftol=posdeftol,drop=FALSE)}))
  sim <- lapply(seq(nrow(sim)), function(i) drop(matrix(sim[i,],
    ncol=ncol(y),byrow=T,dimnames=dimnames(y))))
  if(nsim==1) sim <- sim[[1]]
#
  return(sim)
}
