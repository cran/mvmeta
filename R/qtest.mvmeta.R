qtest.mvmeta <-
function(object, ...) {

  qstat <- object$qstat
	residual <- object$dim$p-object$int>0
	k <- object$dim$k

	res <- list(Q=qstat$Q,df=qstat$df,pvalue=qstat$pvalue,residual=residual,k=k)
	class(res) <- "qtest.mvmeta"
	return(res)
}

