###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012
#
mvmeta <-
  function(formula, S, data, subset, method="reml", model=TRUE, contrasts=NULL,
    offset, na.action, control=list()) {
#
################################################################################
# CREATE THE CALL
#
  call  <- match.call()
  mcall <- match.call(expand.dots=FALSE)
  mn <- match(c("formula","data","constrasts"), names(mcall), 0L)
  mn <- match(c("formula", "data", "subset", "weights", "na.action", 
    "offset"), names(mcall), 0L)
  mcall <- mcall[c(1L, mn)]
  mcall$drop.unused.levels <- TRUE
  mcall[[1L]] <- as.name("model.frame")
#
#  # CREATE FORMULA IF NOT PROVIDED (FOR SIMPLE META-ANALYSIS)
  if(class(eval(mcall[[mn[1]]],
    if(missing(data)) parent.frame() else data))!="formula") {
    formula <- as.formula(paste(deparse(substitute(formula),
      width.cutoff=499L),"~ 1"))
    environment(formula) <- parent.frame()
    call[[mn[1]]] <- mcall[[mn[1]]] <- formula
  }
#
################################################################################
# DERIVE THE MODEL FRAME (SPECIAL HANDLING OF MISSING VALUES)
#
  # NOW KEEP THE MISSING
  mcall$na.action <- "na.pass"
  # CREATE MODEL FRAME WITH ADDITIONAL CLASS
  mf <- eval(mcall,parent.frame())
  class(mf) <- c("data.frame.mvmeta",class(mf))
  # NOW HANDLE THE MISSING
  if(missing(na.action)) na.action <- getOption("na.action")
  if(length(na.action)) mf <- do.call(na.action,list(mf))
  # RETURN mf IF REQUIRED
  method <- match.arg(method,c("fixed","ml","reml","model.frame"))
  if(method=="model.frame") return(mf)
  # EMPTY MODEL?
  if(is.empty.model(mf)) stop("empty model not allowed")
#
################################################################################
# DERIVE OBJECTS FOR FITTING
#
  terms <- attr(mf,"terms")
  # KEEP RESPONSE AS MATRIX
  y <- as.matrix(model.response(mf,"numeric"))
  X <- model.matrix(terms,mf,contrasts)
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y)) 
      stop("number of offsets should equal number of observations")
  }
  # CREATE S (AS A MATRIX OF VECTORIZED COVAR)
  S <- eval(call$S,if(missing(data)) parent.frame() else data)
  S <- .mkS(S,y,attr(mf,"na.action"),if(missing(subset)) NULL else subset)
  if(nrow(y)<2L) stop("less than 2 valid studies after exclusion of missing")
#
################################################################################
# FIT THE MODEL CALLING mvmeta.fit
#  
  # MODEL FIT
  fit <- mvmeta.fit(X,y,S,offset,method,control)
#
################################################################################
#  COMPLETE THE LIST
#
  fit$S <- S
  fit$model <- if(model) mf else NULL
  fit$na.action <- attr(mf,"na.action")
  fit$call <- call
  fit$formula <- formula
  fit$terms <- terms
  fit$offset <- offset 
  fit$control <- control
  fit$method <- method
  fit$contrasts <- attr(X,"contrasts")
  fit$xlevels <- .getXlevels(terms,mf)
#
  class(fit) <- "mvmeta"
#
  return(fit)
}

