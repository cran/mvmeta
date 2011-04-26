meanlist <-
function(list) {
	n <- length(list)
	res <- 0
	for (i in seq(n)) res <- res+list[[i]]
	return(res/n)
}

