xpndMat <-
function(vech) {
	dim <- (-1 + sqrt(1 + 8 * length(vech)))/2
	mat <- matrix(nrow=dim,ncol=dim)
	mat[lower.tri(mat,diag=TRUE)] <- as.matrix(vech)
	mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
	return(mat)
}

