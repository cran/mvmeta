simulate.mvmeta <-
function(beta=rep(1,k), m=100, k=2, Psicomp=c(1,0.2), 
	Scomp=c(1,0.2), Sdf=30) {

	# BETWEEN COVAR MATRIX
	Psi <- matrix(Psicomp[2],k,k)
	diag(Psi) <- diag(Psi)+diff(rev(Psicomp))

	# TRUE CITY-SPECIFIC EFFECTS
	ylist <- lapply(1:m,function(i) {
		as.vector(rmvnorm(n=1,beta,Psi))
	})

	# KNOWN WITHIN COVAR MATRICES (FORCED TO BE POSITIVE)
	S <- matrix(Scomp[2],k,k)
	diag(S) <- diag(S)+diff(rev(Scomp))
	Slist <- lapply(1:m, function(x) rwish(Sdf,S)/Sdf)

	# ESTIMATED CITY-SPECIFIC EFFECTS (WITH IDENTICAL s)
	ylist <- mapply(function(y,S) as.vector(rmvnorm(1,y,S)),ylist,Slist,
		SIMPLIFY=FALSE)

	return(list(beta=beta,S=S,Psi=Psi,ylist=ylist,Slist=Slist))
}

