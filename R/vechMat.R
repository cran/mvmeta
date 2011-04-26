vechMat <-
function(mat, diag=TRUE) {
	if(!is.matrix(mat)) mat <- as.matrix(mat)
	return(mat[lower.tri(mat,diag=diag)])
}

