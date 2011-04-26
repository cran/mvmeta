gradchol <-
function(i,k,U) {
	ind1 <- rep(1:k,k:1)
	ind2 <- unlist(sapply(1:k,seq,to=k))
	C <- D <- E <- diag(0,k)
	C[ind2[i],] <- D[,ind2[i]] <- U[ind1[i],]
	E[ind2[i],] <- E[,ind2[i]] <- 1
	F <- E*C+E*D
	return(F)
}

