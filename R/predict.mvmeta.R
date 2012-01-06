predict.mvmeta <-
function(object, newdata, se=FALSE, ci=FALSE, vcov=FALSE,
	interval=c("confidence","prediction"), ci.level=0.95,
	format=c("matrix","list"), aggregate=c("stat","y"), na.action, ...) {

	# CHECKS
	interval <- match.arg(interval,c("confidence","prediction"))
	if(ci.level<=0||ci.level>=1) stop("'ci.level' must be within 0 and 1")
	format <- match.arg(format,c("matrix","list"))
	aggregate <- match.arg(aggregate,c("stat","y"))
	if(missing(na.action)) na.action <- object$na.action
	na.action <- match.arg(na.action,c("na.omit","na.exclude",
		"na.fail","na.pass"))

	# CREATE DESIGN MATRIX X
	# IF ONLY-INTERCEPT MODEL, ONLY 1 PREDICTED VALUE
	if(object$dim$p-object$int==0) { 
		X <- as.matrix(1)
	# IF newdata MISSING, USE ORIGINAL SET OF PREDICTORS
	} else if(missing(newdata) || is.null(newdata)) {
		model <- object$model
		terms <- terms(model)
		X <- model.matrix(terms,model,object$contrast)
	# IF NOT MISSING
	} else {
		model <- object$model
		terms <- delete.response(terms(model))
		xlevels <- .getXlevels(terms,model)
		# NEW DATA MAY BE ALSO A MATRIX OR VECTOR
		mdata <- model.frame(terms,newdata,
			na.action=na.action,xlev=xlevels)
		X <- model.matrix(terms,mdata,contrasts.arg=object$contrasts)
	}

#########################################################################

	# RE-CREATE OBJECTS
	# HERE nalist INCLUDES ALL THE OBS, ALSO MISSING
	# kXlist IS OF LENGTH 1 IF NO PREDICTOR
	len <- nrow(X)
	nalist <- lapply(seq(len),function(x) rep(TRUE,object$dim$k))
	kXlist <- mapply(function(i,na) {diag(1,object$dim$k)[na,,drop=FALSE]%x%
		X[i,,drop=FALSE]},seq(len),nalist,SIMPLIFY=FALSE)

	# COMPUTE PREDICTION
	predlist <- lapply(kXlist,function(kX) {
		fit <- as.numeric(kX%*%object$coef)
		names(fit) <- object$lab$klab
		return(fit)})

	# COMPUTE VCOV, SE AND CI
	zvalci <- qnorm((1-ci.level)/2,lower.tail=FALSE)
	vcovlist <- lapply(kXlist,function(kX) {
		vcov <- kX%*%tcrossprod(object$vcov,kX)
		if(interval=="prediction") vcov <- vcov+object$Psi
		dimnames(vcov) <- list(object$lab$klab,object$lab$klab)
		return(vcov)})
	selist <- lapply(vcovlist, function(x) sqrt(diag(x)))
	cilblist <- mapply(function(fit,se) fit-zvalci*se,
			predlist,selist,SIMPLIFY=FALSE)
	ciublist <- mapply(function(fit,se) fit+zvalci*se,
			predlist,selist,SIMPLIFY=FALSE)

	# HANDLE MISSING
	Xna <- rowSums(!is.na(X))==object$dim$p
	if(na.action=="na.fail" && !all(Xna)) {
		stop("missing values in 'X'")
	}
	if(na.action=="na.omit") {
		predlist <- predlist[Xna]
		selist <- selist[Xna]
		cilblist <- cilblist[Xna]
		ciublist <- ciublist[Xna]
		vcovlist <- vcovlist[Xna]
	}

#########################################################################

	# IF VCOV WHEN MORE THAN 1 OUTCOME, SWITCH TO LIST
	if(vcov && object$dim$k>1) format <- "list"
	# IF ONLY 1 OUTCOME, FORCE THE AGGREGATION ON IT
	if(format=="matrix" && object$dim$k==1) {
		fit <- rbindlist(predlist)
		if(se) fit <- cbind(fit,rbindlist(selist))
		if(ci) fit <- cbind(fit,rbindlist(cilblist),rbindlist(ciublist))
		if(vcov) fit <- cbind(fit,rbindlist(vcovlist))
		dimnames(fit) <- list(NULL,c("fit","se","ci.lb",
			"ci.ub","vcov")[c(TRUE,se,ci,ci,vcov)])
	# IF NO OTHER STAT, FORCE THE AGGREGATION ON PREDICTION
	}else if(format=="matrix" && !se && !ci) {
		fit <- rbindlist(predlist)
	# WHEN AGGREGATE ON STAT
	} else if (format=="matrix" && aggregate=="stat") {
		fit <- list(fit=rbindlist(predlist))
		if(se) fit$se <- rbindlist(selist)
		if(ci) {
			fit$ci.lb <- rbindlist(cilblist)
			fit$ci.ub <- rbindlist(ciublist)
		}
	# WHEN AGGREGATE ON OUTCOME
	} else if (format=="matrix" && aggregate=="y") {
		fit <- list()
		for(j in seq(object$dim$k)) {
			templist <- lapply(seq(predlist),function(i) {
				fit <- predlist[[i]][j]
				if(se) fit <- cbind(fit,selist[[i]][j])
				if(ci) fit <- cbind(fit,cilblist[[i]][j],
					ciublist[[i]][j])
				return(as.matrix(fit))})
			temp <- rbindlist(templist)
			dimnames(temp) <- list(NULL,
				c("fit","se","ci.lb","ci.ub")[c(TRUE,se,ci,ci)])
			fit[[j]] <- temp
		}
		names(fit) <- object$lab$klab
	# ALL THE OTHER COMBINATIONS, PRODUCE A LIST
	} else {
		fit <- mapply(function(predarg,searg,cilbarg,ciubarg,vcovarg) {
			temp <- list(fit=predarg)
			if(se) temp$se <- searg
			if(ci) temp$ci.lb <- cilbarg
			if(ci) temp$ci.ub <- ciubarg
			if(vcov) temp$vcov <- vcovarg
      # RETURN ONE OBJECT IF ONLY FIT 
			if(se||ci||vcov) {
        return(temp)
      } else return(temp[[1]])},
			predlist,selist,cilblist,ciublist,vcovlist,SIMPLIFY=FALSE)
		# IF PREDICTION ONLY FOR 1 VALUE, ELIMINATE THE FIRST HIERARCHY
		if(length(fit)==1) fit <- fit[[1]]	
	}
  
  # SIMPLIFY IF MATRIX AND ONLY 1 PREDICTION
  if(format=="matrix"  && length(predlist)==1L) {
    if(is.matrix(fit)) fit <- fit[1,]
    if(is.list(fit)) fit <- lapply(fit,function(x) x[1,])
  }

	return(fit)
}

