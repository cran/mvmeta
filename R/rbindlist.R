rbindlist <-
function(list) {
	n <- length(list)
	res <- NULL
	for (i in seq(n)) res <- rbind(res, list[[i]])
	return(res)
}

