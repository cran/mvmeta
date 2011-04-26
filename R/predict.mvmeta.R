predict.mvmeta <-
function(object, newdata, se=FALSE, ci=FALSE, vcov=FALSE,
	interval=c("confidence","prediction"), ci.level=0.95,
	format=c("matrix","list"), aggregate=c("stat","y"), na.action, ...) {

	# CHECKS
	if(is.null(object$X)) {
		if(!missing(newdata)) stop("original model does not include predictors")
		newdata <- NULL
	}
	if(missing(newdata)) newdata <- object$X
	if(is.list(newdata)&&!is.data.frame(newdata)) newdata <- rbindlist(newdata)
	if(!is.null(newdata) && is.data.frame(newdata)) {
		newdata <- as.matrix(newdata)
	}
	if(!is.null(newdata) && !is.matrix(newdata)) {
		if(object$dim$p==2) { newdata <- as.matrix(newdata)
		} else newdata <- t(as.matrix(newdata))
	}
	if(!is.null(newdata) && ncol(newdata)!=object$dim$p-1) {
		stop("columns in 'newdata' must match number of predictors")
	}
	interval <- match.arg(interval,c("confidence","prediction"))
	if(ci.level<=0||ci.level>=1) stop("'ci.level' must be within 0 and 1")
	format <- match.arg(format,c("matrix","list"))
	aggregate <- match.arg(aggregate,c("stat","y"))
	if(missing(na.action)) na.action <- object$na.action
	na.action <- match.arg(na.action,c("na.omit","na.exclude",
		"na.fail","na.pass"))

#########################################################################

	# RE-CREATE OBJECTS
	# HERE nalist INCLUDES ALL THE OBS, ALSO MISSING
	# kXlist IS OF LENGTH 1 IF NO PREDICTOR
	len <- ifelse(!is.null(newdata),nrow(newdata),1)
	nalist <- lapply(seq(len),function(x) rep(TRUE,object$dim$k))
	kXlist <- with(object,kXlistmk(newdata,cen,nalist,len,dim$k))

	# COMPUTE PREDICTION
	predlist <- lapply(kXlist,function(kX) {
		pred <- as.numeric(kX%*%object$beta)
		names(pred) <- object$lab$klab
		return(pred)})

	# COMPUTE VCOV, SE AND CI
	zvalci <- qnorm((1-ci.level)/2,lower.tail=FALSE)
	vcovlist <- lapply(kXlist,function(kX) {
		vcov <- kX%*%tcrossprod(object$vcov,kX)
		if(interval=="prediction") vcov <- vcov+object$Psi
		dimnames(vcov) <- list(object$lab$klab,object$lab$klab)
		return(vcov)})
	selist <- lapply(vcovlist, function(x) sqrt(diag(x)))
	cilblist <- mapply(function(pred,se) pred-zvalci*se,
			predlist,selist,SIMPLIFY=FALSE)
	ciublist <- mapply(function(pred,se) pred+zvalci*se,
			predlist,selist,SIMPLIFY=FALSE)

	# HANDLE MISSING
	if(!is.null(newdata)) {
		Xna <- rowSums(!is.na(object$X))==ncol(object$X)
		if(na.action=="na.fail" && !all(Xna)) {
			stop("missing values in 'X'")
		}
		if(!all(Xna) && na.action=="na.omit") {
			predlist <- predlist[Xna]
			selist <- selist[Xna]
			cilblist <- cilblist[Xna]
			ciublist <- ciublist[Xna]
			vcovlist <- vcovlist[Xna]
		}
	}

#########################################################################

	# IF VCOV WHEN MORE THAN 1 OUTCOME, SWITCH TO LIST
	if(vcov && object$dim$k>1) format <- "list"
	# IF ONLY 1 OUTCOME, FORCE THE AGGREGATION ON IT
	if(format=="matrix" && object$dim$k==1) {
		pred <- rbindlist(predlist)
		if(se) pred <- cbind(pred,rbindlist(selist))
		if(ci) pred <- cbind(pred,rbindlist(cilblist),rbindlist(ciublist))
		if(vcov) pred <- cbind(pred,rbindlist(vcovlist))
		dimnames(pred) <- list(NULL,c("pred","se","ci.lb",
			"ci.ub","vcov")[c(TRUE,se,ci,ci,vcov)])
	# IF NO OTHER STAT, FORCE THE AGGREGATION ON PREDICTION
	}else if(format=="matrix" && !se && !ci) {
		pred <- rbindlist(predlist)
	# WHEN AGGREGATE ON STAT
	} else if (format=="matrix" && aggregate=="stat") {
		pred <- list(pred=rbindlist(predlist))
		if(se) pred$se <- rbindlist(selist)
		if(ci) {
			pred$ci.lb <- rbindlist(cilblist)
			pred$ci.ub <- rbindlist(ciublist)
		}
	# WHEN AGGREGATE ON OUTCOME
	} else if (format=="matrix" && aggregate=="y") {
		pred <- list()
		for(j in seq(object$dim$k)) {
			templist <- lapply(seq(predlist),function(i) {
				pred <- predlist[[i]][j]
				if(se) pred <- cbind(pred,selist[[i]][j])
				if(ci) pred <- cbind(pred,cilblist[[i]][j],
					ciublist[[i]][j])
				return(as.matrix(pred))})
			temp <- rbindlist(templist)
			dimnames(temp) <- list(NULL,
				c("pred","se","ci.lb","ci.ub")[c(TRUE,se,ci,ci)])
			pred[[j]] <- temp
		}
		names(pred) <- object$lab$klab
	# ALL THE OTHER COMBINATIONS, PRODUCE A LIST
	} else {
		pred <- mapply(function(predarg,searg,cilbarg,ciubarg,vcovarg) {
			temp <- list(pred=predarg)
			if(se) temp$se <- searg
			if(ci) temp$ci.lb <- cilbarg
			if(ci) temp$ci.ub <- ciubarg
			if(vcov) temp$vcov <- vcovarg
			return(temp)},
			predlist,selist,cilblist,ciublist,vcovlist,SIMPLIFY=FALSE)
		# IF PREDICTION ONLY FOR 1 VALUE, ELIMINATE THE FIRST HIERARCHY
		if(length(pred)==1) pred <- pred[[1]]	
	}

	return(pred)
}

