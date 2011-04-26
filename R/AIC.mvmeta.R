AIC.mvmeta <-
function (object, ..., k=2) {
	AIC(logLik(object), k=k)
}

