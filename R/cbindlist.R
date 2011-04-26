cbindlist <-
function(list) {
	n <- length(list)
	res <- NULL
	for (i in seq(n)) res <- cbind(res, list[[i]])
	return(res)
}

