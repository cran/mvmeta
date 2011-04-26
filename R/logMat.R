logMat <-
function(x, positive=FALSE) {
	if(is.data.frame(x)) x <- as.matrix(x)
	if(is.matrix(x)) {
		expA <- x
	} else expA <- xpndMat(x)
	eig <- eigen(expA)
	if(positive) eig$val[eig$val<10^-9] <- 10^-9
	B <- diag(log(eig$val),nrow(expA))
	A <- eig$vec%*%tcrossprod(B,eig$vec)
	return(A)
}

