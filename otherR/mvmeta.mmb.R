###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2013
#
`mvmeta.mmb` <-
function(Xlist, ylist, Slist, nalist, k, m, p, nall, control) {
#
################################################################################
#
  # FIT FIXED EFFECTS MODEL
  Psi <- diag(0,k)
  gls <- .gls(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#
  Qlist <- mapply(function(invU,y,X) tcrossprod(invU) %*% 
    tcrossprod(y-X%*%gls$coef),gls$invUlist,ylist,Xlist,SIMPLIFY=FALSE)
  Q <- matrix(0,nrow=k,ncol=k)
  for(i in seq(m)) Q[!nalist[[i]],!nalist[[i]]] <- 
    Q[!nalist[[i]],!nalist[[i]]] + Qlist[[i]]
#
  # COMPUTE A AND B LISTS WITH BLOCKS
  # NB: Hij <- Xi %*% tXWXtot %*% t(Xj) %*% Wjj
  #     Aij <- t(I-Hji) %*% Wjj 
  Alist <- Blist <- vector("list",m^2)
  tXWXtot <- .sumlist(lapply(gls$invtUXlist,crossprod))
  invtXWXtot <- chol2inv(chol(tXWXtot))
  ind1 <- rep(seq(m),m)
  ind2 <- rep(seq(m),each=m)
  for(i in seq(Alist)) {
    Hji <- Xlist[[ind2[i]]] %*% invtXWXtot %*% t(Xlist[[ind1[i]]]) %*%
      tcrossprod(gls$invUlist[[ind1[i]]])
    diag <- diag(ind1[[i]]==ind2[[i]],k)
    diag <- diag[!nalist[[ind2[i]]],!nalist[[ind1[i]]],drop=FALSE]
    Alist[[i]] <- crossprod(diag-Hji,tcrossprod(gls$invUlist[[ind2[i]]]))
    Blist[[i]] <- t(diag-Hji)
  }
#
  # COMPUTE COMPONENTS OF SYSTEM OF EQUATIONS
  Bcomp <- matrix(0,nrow=k,ncol=k)
  for(i in diag(matrix(seq(m*m),nrow=m))) {
    Bcomp[!nalist[[ind1[i]]],!nalist[[ind2[i]]]] <- 
    Bcomp[!nalist[[ind1[i]]],!nalist[[ind2[i]]]] + Blist[[i]]
  }
  # NB: Aij = t(Aji), SO t(Bji)%x%Aij = t(Bji%x%Aji)
  tBkAcomp <- matrix(0,nrow=k*k,ncol=k*k)
  for(i in seq(Alist)) {
    mat <- matrix(0,k,k)
    mat[!nalist[[ind2[i]]],!nalist[[ind1[i]]]] <- 1
    mat <- mat%x%mat
    tBkAcomp[mat==1] <- tBkAcomp[mat==1] + t(Blist[[i]] %x% Alist[[i]])
  }   
  # SOLVE THE SYSTEM
  Psi1 <- qr.solve(tBkAcomp,as.numeric(Q-Bcomp))
#
  # CREATE Psi
  dim(Psi1) <- c(k,k)
  Psi <- (Psi1+t(Psi1))/2
#
  # FORCE SEMI-POSITIVE DEFINITENESS
  eig <- eigen(Psi)
  Psi <- eig$vectors %*% diag(pmax(eig$values,0),k) %*% t(eig$vectors)
#
  # FIT BY GLS
  gls <- .gls(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#
  # COMPUTE (CO)VARIANCE MATRIX OF coef
  qrinvtUX <- qr(gls$invtUX)
  R <- qr.R(qrinvtUX)
  Qty <- qr.qty(qrinvtUX,gls$invtUy)
  vcov <- tcrossprod(backsolve(R,diag(1,ncol(gls$invtUX))))
#
  # COMPUTE RESIDUALS (LATER), FITTED AND RANK
  res <- NULL
  fitted <- lapply(Xlist,"%*%",gls$coef)
  rank <- qrinvtUX$rank
#
  fit <- list(coefficients=gls$coef,vcov=vcov,Psi=Psi,residuals=res,
    fitted.values=fitted,df.residual=nall-rank-length(par),rank=rank,logLik=NA,
    control=control)
#
  return(fit)
}
