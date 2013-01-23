###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`vechMat` <-
function(mat, diag=TRUE) {
#
################################################################################
#
  if(!is.matrix(mat)) mat <- as.matrix(mat)
  if (dim(mat)[1]!=dim(mat)[2]) stop("Non-square matrix")
  return(mat[lower.tri(mat,diag=diag)])
}

