expMat <-
function(x) {
	if(is.data.frame(x)) x <- as.matrix(x)
	if(is.matrix(x)) {
		logA <- x
	} else logA <- xpndMat(x)
	eig <- eigen(logA)
	B <- diag(exp(eig$val),nrow(logA))
	A <- eig$vec%*%tcrossprod(B,eig$vec)
	return(A)
}

