
#SCCS 12/22/95 @(#)Surv.s	4.17
# Package up surivival type data as a structure
#
Surv <- function(time, time2, event,
	      type=c('right', 'left', 'interval', 'counting', 'interval2'),
		       origin=0) {
    nn <- length(time)
    ng <- nargs()
    if (missing(type)) {
	if (ng==1 || ng==2) type <- 'right'
	else if (ng==3)     type <- 'counting'
	else stop("Invalid number of arguments")
	}
    else {
	type <- match.arg(type)
	ng <- ng-1
	if (ng!=3 && (type=='interval' || type =='counting'))
		stop("Wrong number of args for this type of survival data")
	if (ng!=2 && (type=='right' || type=='left' ||  type=='interval2'))
		stop("Wrong number of args for this type of survival data")
	}
    who <- !is.na(time)

    if (ng==1) {
	if (!is.numeric(time)) stop ("Time variable is not numeric")
	else if (any(time[who]<0))  stop ("Time variable must be >= 0")
	ss <- cbind(time, 1)
	dimnames(ss) <- list(NULL, c("time", "status"))
	}
    else if (type=='right' || type=='left') {
	if (!is.numeric(time)) stop ("Time variable is not numeric")
	else if (any(time[who]<0))  stop ("Time variable must be >= 0")
	if (length(time2) != nn) stop ("Time and status are different lengths")
	if (is.logical(time2)) status <- 1*time2
	    else  if (is.numeric(time2)) {
		who2 <- !is.na(time2)
                ###<TSL> don't want all 0 to come up as all failed.
                if (max(time2[who2])>1) 
		    status <- time2 - (max(time2[who2]) -1)
                else
                    status<-time2
                ###<TSL>
		if (any(status[who2] !=0  & status[who2]!=1))
				stop ("Invalid status value")
		}
	    else stop("Invalid status value")
	 ss <- cbind(time, status)
	 dimnames(ss) <- list(NULL, c("time", "status"))
	}
    else  if (type=='counting') {
	if (length(time2) !=nn) stop ("Start and stop are different lengths")
	if (length(event)!=nn) stop ("Start and event are different lengths")
	if (!is.numeric(time))stop("Start time is not numeric")
	if (!is.numeric(time2)) stop("Stop time is not numeric")
	who3 <- who & !is.na(time2)
	if (any (time[who3]>= time2[who3]))stop("Stop time must be > start time")
	if (is.logical(event)) status <- 1*event
	    else  if (is.numeric(event)) {
		who2 <- !is.na(event)
		status <- event - min(event[who2])
		if (any(status[who2] !=0  & status[who2]!=1))
				stop("Invalid status value")
		}
	    else stop("Invalid status value")
	ss <- cbind(time-origin, time2-origin, status)
	}

    else {  #interval censored data
	if (type=='interval2') {
	    event <- ifelse(is.na(time), 2,
		     ifelse(is.na(time2),0,
		     ifelse(time==time2, 1,3)))
	    if (any(time[event==3] > time2[event==3]))
		stop("Invalid interval: start > stop")
	    time <- ifelse(event!=2, time, time2)
	    type <- 'interval'
	    }
	else {
	    temp <- event[!is.na(event)]
	    if (!is.numeric(temp)) stop("Status indicator must be numeric")
	    if (length(temp)>0 && any(temp!= floor(temp) | temp<0 | temp>3))
		stop("Status indicator must be 0, 1, 2 or 3")
	    }
	status <- event
	ss <- cbind(time, ifelse(!is.na(event) & event==3, time2, 1),
			    status)
	}

    attr(ss, "class") <- c("Surv")
    attr(ss, "type")  <- type
    ss
    }

print.Surv <- function(xx, quote=F, ...)
    invisible(print(as.character.Surv(xx), quote=quote, ...))

as.character.Surv <- function(xx) {
    class(xx) <- NULL
    type <- attr(xx, 'type')
    if (type=='right') {
	temp <- xx[,2]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+"," "))
	paste(format(xx[,1]), temp, sep='')
	}
    else if (type=='counting') {
	temp <- xx[,3]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+"," "))
	paste('(', format(xx[,1]), ',', format(xx[,2]), temp,
			 ']', sep='')
	}
    else if (type=='left') {
	temp <- xx[,2]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "<"," "))
	paste(temp, format(xx[,1]), sep='')
	}
    else {   #interval type
	stat <- xx[,3]
	temp <- c("+", "", "-", "]")[stat+1]
	temp2 <- ifelse(stat==3,
			 paste("[", format(xx[,1]), ", ",format(xx[,2]), sep=''),
			 format(xx[,1]))
	ifelse(is.na(stat), "NA", paste(temp2, temp, sep=''))
	}
    }

"[.Surv" <- function(x,i,j, drop=F) {
   temp <- class(x)
    type <- attr(x, "type")
    class(x) <- NULL
    attr(x, 'type') <- NULL
    if (missing(j)) {
	x <- x[i,,drop=drop]
	class(x) <- temp
	attr(x, "type") <- type
	x
	}
    else NextMethod("[")
    }

is.na.Surv <- function(x) {
    class(x) <- NULL
    as.vector( (1* is.na(x))%*% rep(1, ncol(x)) >0)
    }

Math.Surv <- function(...)  stop("Invalid operation on a survival time")
Ops.Surv  <- function(...)  stop("Invalid operation on a survival time")
Summary.Surv<- function(...) stop("Invalid operation on a survival time")
is.Surv <- function(x) inherits(x, 'Surv')
#SCCS @(#)agexact.fit.s	4.16 02/28/95
agexact.fit <- function(x, y, strata, offset, iter.max,
			eps, weights, init, method, rownames)
    {
    if (!is.matrix(x)) stop("Invalid formula for cox fitting function")
    if (!is.null(weights) && any(weights!=1))
	  stop("Case weights are not supported for the exact method")
    n <- nrow(x)
    nvar <- ncol(x)
    if (ncol(y)==3) {
	start <- y[,1]
	stopp <- y[,2]
	event <- y[,3]
	}
    else {
	start <- rep(0,n)
	stopp <- y[,1]
	event <- y[,2]
	}

    # Sort the data (or rather, get a list of sorted indices)
    if (length(strata)==0) {
	sorted <- order(stopp, -event)
	newstrat <- as.integer(rep(0,n))
	}
    else {
	sorted <- order(strata, stopp, -event)
        ## change as.numeric() to codes()
	strata <-(as.numeric(strata))[sorted]
	newstrat <- as.integer(c(1*(diff(strata)!=0), 1))
	}
    if (is.null(offset)) offset <- rep(0,n)

    sstart <- as.double(start[sorted])
    sstop <- as.double(stopp[sorted])
    sstat <- as.integer(event[sorted])

    if (is.null(nvar)) {
	# A special case: Null model.  Not worth coding up
	stop("Cannot handle a null model + exact calculation (yet)")
	}

    if (!is.null(init)) {
	if (length(init) != nvar) stop("Wrong length for inital values")
	}
    else init <- rep(0,nvar)

    agfit <- .C("agexact", iter= as.integer(iter.max),
		   as.integer(n),
		   as.integer(nvar), sstart, sstop,
		   sstat,
		   x= x[sorted,],
		   as.double(offset[sorted] - mean(offset)),
		   newstrat,
		   means = double(nvar),
		   coef= as.double(init),
		   u = double(nvar),
		   imat= double(nvar*nvar), loglik=double(2),
		   flag=integer(1),
		   double(2*nvar*nvar +nvar*4 + n),
		   integer(2*n),
		   as.double(eps),
		   sctest=double(1) )

    var <- matrix(agfit$imat,nvar,nvar)
    coef <- agfit$coef
    if (agfit$flag < nvar) which.sing <- diag(var)==0
    else which.sing <- rep(F,nvar)

    infs <- abs(agfit$u %*% var)
    if (iter.max >1) {
	if (agfit$flag == 1000)
	       warning("Ran out of iterations and did not converge")
	else if (any((infs > eps) & (infs > sqrt(eps)*abs(coef))))
	    warning(paste("Loglik converged before variable ",
			  paste((1:nvar)[(infs>eps)],collapse=","),
			  "; beta may be infinite. "))
	}

    names(coef) <- dimnames(x)[[2]]
    lp  <- x %*% coef + offset - sum(coef *agfit$means)
    score <- as.double(exp(lp[sorted]))
    agres <- .C("agmart",
		   as.integer(n),
		   as.integer(0),
		   sstart, sstop,
		   sstat,
		   score,
		   rep(1,n),
		   newstrat,
		   resid=double(n))
    resid _ double(n)
    resid[sorted] <- agres$resid
    names(resid) <- rownames
    coef[which.sing] <- NA

    list(coefficients  = coef,
		var    = var,
		loglik = agfit$loglik,
		score  = agfit$sctest,
		iter   = agfit$iter,
		linear.predictors = lp,
		residuals = resid,
		means = agfit$means,
		method= 'coxph')
    }
#SCCS @(#)agreg.fit.s	4.15 02/28/95
agreg.fit <- function(x, y, strata, offset, init, iter.max,
			eps, weights, method, rownames)
    {
    n <- nrow(y)
    nvar <- ncol(x)
    start <- y[,1]
    stopp <- y[,2]
    event <- y[,3]

    # Sort the data (or rather, get a list of sorted indices)
    if (length(strata)==0) {
	sorted <- order(stopp, -event)
	newstrat <- as.integer(rep(0,n))
	}
    else {
	sorted <- order(strata, stopp, -event)
        ## change as.numeric() to codes()
	strata <- (as.numeric(strata))[sorted]
	newstrat <- as.integer(c(1*(diff(strata)!=0), 1))
	}
    if (missing(offset) || is.null(offset)) offset <- rep(0,n)
    if (missing(weights)|| is.null(weights))weights<- rep(1,n)
    else {
	if (any(weights<=0)) stop("Invalid weights, must be >0")
	weights <- weights[sorted]
	}
    sstart <- as.double(start[sorted])
    sstop <- as.double(stopp[sorted])
    sstat <- as.integer(event[sorted])

    if (is.null(nvar)) {
	# A special case: Null model.  Just return obvious stuff
	score <- as.double(exp(offset[sorted]))
	agfit <- .C("agfit_null",
		       as.integer(n),
		       as.integer(method=='efron'),
		       sstart, sstop,
		       sstat,
		       offset[sorted],
		       as.double(weights),
		       newstrat,
		       loglik=double(1))

	agres <- .C("agmart",
		       as.integer(n),
		       as.integer(method=='efron'),
		       sstart, sstop,
		       sstat,
		       score,
		       as.double(weights),
		       newstrat,
		       resid=double(n))

	resid _ double(n)
	resid[sorted] <- agres$resid
	names(resid) <- rownames

	list(loglik=agfit$loglik,
	     linear.predictors = offset,
	     residuals = resid,
	     method= c("coxph.null", 'coxph') )
	}

    else {
	if (!is.null(init)) {
	    if (length(init) != nvar) stop("Wrong length for inital values")
	    }
	else init <- rep(0,nvar)

	agfit <- .C("agfit2", iter= as.integer(iter.max),
		       as.integer(n),
		       as.integer(nvar), sstart, sstop,
		       sstat,
		       x= x[sorted,],
		       as.double(offset[sorted] - mean(offset)),
		       as.double(weights),
		       newstrat,
		       means = double(nvar),
		       coef= as.double(init),
		       u = double(nvar),
		       imat= double(nvar*nvar), loglik=double(2),
		       flag=integer(1),
		       double(2*nvar*nvar +nvar*3 + n),
		       integer(n),
		       as.double(eps),
		       sctest=as.double(method=='efron') )

	var <- matrix(agfit$imat,nvar,nvar)
	coef <- agfit$coef
	if (agfit$flag < nvar) which.sing <- diag(var)==0
	else which.sing <- rep(F,nvar)

	infs <- abs(agfit$u %*% var)
	if (iter.max >1) {
	    if (agfit$flag == 1000)
		   warning("Ran out of iterations and did not converge")
	    else if (any((infs > eps) & (infs > sqrt(eps)*abs(coef))))
		warning(paste("Loglik converged before variable ",
			  paste((1:nvar)[(infs>eps)],collapse=","),
			  "; beta may be infinite. "))
	    }

	names(coef) <- dimnames(x)[[2]]
	lp  <- x %*% coef + offset - sum(coef *agfit$means)
	score <- as.double(exp(lp[sorted]))
	agres <- .C("agmart",
		       as.integer(n),
		       as.integer(method=='efron'),
		       sstart, sstop,
		       sstat,
		       score,
		       as.double(weights),
		       newstrat,
		       resid=double(n))

	resid _ double(n)
	resid[sorted] <- agres$resid
	names(resid) <- rownames
	coef[which.sing] <- NA

	list(coefficients  = coef,
		    var    = var,
		    loglik = agfit$loglik,
		    score  = agfit$sctest,
		    iter   = agfit$iter,
		    linear.predictors = as.vector(lp),
		    residuals = resid,
		    means = agfit$means,
		    method= 'coxph')
	}
    }
#SCCS @(#)cluster.s	4.1 10/01/94
cluster <- function(x) x
# SCCS @(#)cox.zph.s	1.13 12/27/95
#  Test proportional hazards
#
cox.zph <- function(fit, transform='km', global=T) {
    call <- match.call()
    if (!inherits(fit, 'coxph')) stop ("Argument must be the result of coxph")
    if (inherits(fit, 'coxph.null'))
	stop("The are no score residuals for a Null model")

    sresid <- resid(fit, 'schoenfeld')
    varnames <- names(fit$coef)
    nvar <- length(varnames)
    ndead<- length(sresid)/nvar
    if (nvar==1) times <- as.numeric(names(sresid))
    else         times <- as.numeric(dimnames(sresid)[[1]])

    if (missing(transform) && attr(fit$y, 'type') != 'right')
	    transform <- 'identity'
    if (is.character(transform)) {
	tname <- transform
	ttimes <- switch(transform,
			   'identity'= times,
			   'rank'    = rank(times),
			   'log'     = log(times),
			   'km' = {
				temp <- survfit.km(factor(rep(1,nrow(fit$y))),
						    fit$y, se.fit=F)
				# A nuisance to do left cont KM
				t1 <- temp$surv[temp$n.event>0]
				t2 <- temp$n.event[temp$n.event>0]
				km <- rep(c(1,t1), c(t2,0))
				if (is.null(attr(sresid, 'strata')))
				    1-km
				else (1- km[sort.list(sort.list(times))])
				},
			   stop("Unrecognized transform"))
	}
    else {
	tname <- deparse(substitute(transform))
	ttimes <- transform(times)
	}
    xx <- ttimes - mean(ttimes)
    ### <TSL> in R can't do 1x1 matrix times vector
    if (nvar==1) sresid<-as.matrix(sresid)
    r2 <- sresid %*% fit$var * ndead
    test <- xx %*% r2        # time weighted col sums
    corel <- c(cor(xx, r2))
    z <- c(test^2 /(diag(fit$var)*ndead* sum(xx^2)))
    Z.ph <- cbind(corel, z, 1- pchisq(z,1))

    if (global && nvar>1) {
	test <- c(xx %*% sresid)
	z    <- c(test %*% fit$var %*% test) * ndead / sum(xx^2)
	Z.ph <- rbind(Z.ph, c(NA, z, 1-pchisq(z, ncol(sresid))))
	dimnames(Z.ph) <- list(c(varnames, "GLOBAL"), c("rho", "chisq", "p"))
	}
    else dimnames(Z.ph) <- list(varnames, c("rho", "chisq", "p"))

    dimnames(r2) <- list(times, names(fit$coef))
    temp <-list(table=Z.ph, x=ttimes, y=r2 + outer(rep(1,ndead), fit$coef),
    var=fit$var, call=call, transform=tname)
    class(temp) <- "cox.zph"
    temp
    }

"[.cox.zph" <- function(x,i, ..., drop=F) {
    z<- list(table=x$table[i,,drop=drop], x=x$x, y=x$y[ ,i,drop=drop],
		var=x$var[i,i, drop=drop], call=x$call,
		transform=x$transform)
    attributes(z) <- attributes(x)
    z
    }
#SCCS  @(#)coxph.detail.s	4.10 07/28/95
coxph.detail <-  function(object) {
    method <- object$method
    if (method!='breslow' && method!='efron')
	stop(paste("Detailed output is not available for the", method,
			"method"))
    n <- length(object$residuals)
    rr <- object$residual
    weights <- object$weights        #always present if there are weights
    x <- object$x
    y <- object$y
    strat <- object$strata
    Terms <- object$terms
    if (!inherits(Terms, 'terms'))
	    stop("invalid terms component of object")
    strats <- attr(Terms, "specials")$strata

    if (is.null(y)  ||  is.null(x)) {
	temp <- coxph.getdata(object, y=T, x=T, strata=T)
	y <- temp$y
	x <- temp$x
	if (length(strats)) strat <- temp$strata
	}

    nvar <- ncol(x)
    if (ncol(y)==2) y <- cbind(-1,y)
    if (is.null(strat)) {
	ord <- order(y[,2], -y[,3])
	newstrat <- rep(0,n)
	}
    else {
	ord <- order(strat, y[,2], -y[,3])
        ## change as.numeric() to codes()
	newstrat <- c(diff(as.numeric(strat[ord]))!=0 ,1)
	}
    newstrat[n] <- 1

    # sort the data
    x <- x[ord,]
    y <- y[ord,]
    storage.mode(y) <- 'double'
    score <- exp(object$linear.predictor)[ord]
    if (is.null(weights)) weights <- rep(1,n)
    else                  weights <- weights[ord]

    ndeath <- sum(y[,3])
    ff <- .C("coxdetail", as.integer(n),
			  as.integer(nvar),
			  ndeath= as.integer(ndeath),
			  y = y,
			  as.double(x),
			  as.integer(newstrat),
			  index =as.double(score),
			  weights = as.double(weights),
			  means= c(method=='efron', double(ndeath*nvar)),
			  u = double(ndeath*nvar),
			  i = double(ndeath*nvar*nvar),
			  double(nvar*(3 + 2*nvar)) )
    keep <- 1:ff$ndeath
    vname<- dimnames(x)[[2]]
    time <- y[ff$index[keep],2]
    names(time) <- NULL
    means<- (matrix(ff$means,ndeath, nvar))[keep,]
    score<-  matrix(ff$u, ndeath, nvar)[keep,]
    var <- array(ff$i, c(nvar, nvar, ndeath))[,,keep]
    if (nvar>1) {
	dimnames(means) <- list(time, vname)
	dimnames(score) <- list(time, vname)
	dimnames(var) <- list(vname, vname, time)
	}
    else {
	names(means) <- time
	names(score) <- time
	names(var) <- time
	}

    dimnames(ff$y) <- NULL
    temp <- list(time = time, means=means, nevent=ff$y[keep,1],
	 nrisk = ff$y[keep,2], hazard= ff$y[keep,3], score= score,  imat=var,
	 varhaz=ff$weights[keep], y=y, x=x)
    if (length(strats)) temp$strata <- table((strat[ord])[ff$index[keep]])
    if (!all(weights==1)) temp$weights <- weights
    temp
    }
#SCCS 04/09/96 @(#)coxph.fit.s	4.18
coxph.fit <- function(x, y, strata, offset, init, iter.max,
			eps, weights, method, rownames)
    {
    n <-  nrow(y)
    if (is.matrix(x)) nvar <- ncol(x)
    else  if (length(x)==0) nvar <-0 else nvar <-1
    time <- y[,1]
    status <- y[,2]

    # Sort the data (or rather, get a list of sorted indices)
    if (length(strata)==0) {
	sorted <- order(time)
	newstrat <- as.integer(rep(0,n))
	}
    else {
	sorted <- order(strata, time)
        ## as.numeric() to codes()
	strata <- (as.numeric(strata))[sorted]
	newstrat <- as.integer(c(1*(diff(strata)!=0), 1))
	}
    if (missing(offset) || is.null(offset)) offset <- rep(0,n)
    if (missing(weights)|| is.null(weights))weights<- rep(1,n)
    else {
	if (any(weights<=0)) stop("Invalid weights, must be >0")
	weights <- weights[sorted]
	}
    stime <- as.double(time[sorted])
    sstat <- as.integer(status[sorted])

    if (nvar==0) {
	# A special case: Null model.
	#  (This is why I need the rownames arg- can't use x' names)
	score <- exp(offset[sorted])
	coxfit <- .C("coxfit_null", as.integer(n),
				    as.integer(method=='efron'),
				    stime,
				    sstat,
				    exp(offset[sorted]),
				    as.double(weights),
				    newstrat,
				    loglik=double(1),
				    resid = double(n) )
	resid <- double(n)
	resid[sorted] <- coxfit$resid
	names(resid) <- rownames

	list( loglik = coxfit$loglik,
	      linear.predictors = offset,
	      residuals = resid,
	      method= c('coxph.null', 'coxph') )
	}

    else {
	if (!missing(init) && !is.null(init)) {
	    if (length(init) != nvar) stop("Wrong length for inital values")
	    }
	else init <- rep(0,nvar)

	coxfit <- .C("coxfit2", iter=as.integer(iter.max),
		       as.integer(n),
		       as.integer(nvar), stime,
		       sstat,
		       x= x[sorted,] ,
		       as.double(offset[sorted] - mean(offset)),
		       as.double(weights),
		       newstrat,
		       means= double(nvar),
		       coef= as.double(init),
		       u = double(nvar),
		       imat= double(nvar*nvar), loglik=double(2),
		       flag=integer(1),
		       double(2*n + 2*nvar*nvar + 3*nvar),
		       as.double(eps),
		       sctest=as.double(method=="efron") )

	var <- matrix(coxfit$imat,nvar,nvar)
	coef <- coxfit$coef
	if (coxfit$flag < nvar) which.sing <- diag(var)==0
	else which.sing <- rep(F,nvar)

	infs <- abs(coxfit$u %*% var)
	if (iter.max >1) {
	    if (coxfit$flag == 1000)
		   warning("Ran out of iterations and did not converge")
	    else if (any((infs > eps) & (infs > sqrt(eps)*abs(coef))))
		warning(paste("Loglik converged before variable ",
			  paste((1:nvar)[(infs>eps)],collapse=" ,"),
			  "; beta may be infinite. ", sep=''))
	    }

	names(coef) <- dimnames(x)[[2]]
	lp <- c(x %*% coef) + offset - sum(coef*coxfit$means)
	score <- exp(lp[sorted])
	coxres <- .C("coxmart", as.integer(n),
				as.integer(method=='efron'),
				stime,
				sstat,
				newstrat,
				score,
				weights,
				resid=double(n))
	resid <- double(n)
	resid[sorted] <- coxres$resid
	names(resid) <- rownames
	coef[which.sing] <- NA

	list(coefficients  = coef,
		    var    = var,
		    loglik = coxfit$loglik,
		    score  = coxfit$sctest,
		    iter   = coxfit$iter,
		    linear.predictors = as.vector(lp),
		    residuals = resid,
		    means = coxfit$means,
		    method='coxph')
	}
    }
# SCCS @(#)coxph.getdata.s	4.3 09/08/94
#
# Reconstruct the Cox model data.  This is done in so many routines
#  that I extracted it out.
#
# The "stratax" name is to avoid conflicts with the strata() function, but
#   still allow users to type "strata" as an arg.
#
coxph.getdata <- function(fit, y=T, x=T, stratax=T, offset=F) {
    ty <- fit$y
    tx <- fit$x
    strat <- fit$strata
    Terms <- fit$terms
    if (is.null(attr(Terms, 'offset'))) offset <- F
    if (offset) x<- T
    if (!inherits(Terms, 'terms'))
	    stop("invalid terms component of fit")
    strats <- attr(Terms, "specials")$strata
    if (length(strats)==0) stratax <- F

    if ( (y && is.null(ty)) || (x && is.null(tx)) ||
	     (stratax && is.null(strat)) || offset) {
	# get the model frame
	m <- fit$model
	if (is.null(m)) m <- model.frame(fit)

	# Pull things out
	if (y && is.null(ty)) ty <- model.extract(m, 'response')

	if (offset) toff <- model.extract(m, 'offset')

	# strata was saved in the fit if and only if x was
	if (x && is.null(tx)) {
	    dropx <- untangle.specials(Terms, 'cluster')$terms
	    if (stratax) {
		temp <- untangle.specials(Terms, 'strata', 1)
		tx <- model.matrix("[.terms"(Terms,-c(dropx,temp$terms)), m)[,-1,drop=F]
		strat <- strata(m[temp$vars], shortlabel=T)
		}
	    else {
		if (length(dropx)) tx <- model.matrix("[.terms"(Terms,-dropx), m)[,-1,drop=F]
		else               tx <- model.matrix(Terms, m)[,-1,drop=F]
		}
	    }
	}
    else if (offset)
       toff <- fit$linear.predictors -(c(tx %*% fit$coef) - sum(fit$means*fit$coef))

    temp <- NULL
    if (y) temp <- c(temp, "y=ty")
    if (x) temp <- c(temp, "x=tx")
    if (stratax)  temp <- c(temp, "strata=strat")
    if (offset)  temp <- c(temp, "offset=toff")

###    eval(parse(text=paste("list(", paste(temp, collapse=','), ")")))
    eval(parse(text=paste("list(", paste(temp, collapse=','), ")"))[[1]])
    }
#SCCS @(#)coxph.rvar.s	4.3 01/18/94
coxph.rvar <- function(fit, collapse) {
    rcall <- match.call()
    if (class(fit) != 'coxph')
	stop ("First argument must be a fitted Cox model")

    if (missing(collapse)) temp <- residuals.coxph(fit, type='dfbeta')
    else temp <- residuals.coxph(fit, type='dfbeta', collapse=collapse)
    if (any(is.na(temp)))
       if (ncol(temp)==1) temp<- temp[!is.na(temp),,drop=F]
       else               temp <- temp[!is.na(temp %*% rep(1,ncol(temp))),]
    fit$robust.var <- t(temp) %*% temp
    fit$rcall <- rcall
    fit
    }
#SCCS  @(#)coxph.s	4.20  08/30/96
coxph <- function(formula=formula(data), data=sys.frame(sys.parent()),
	weights, subset, na.action,
	eps=.0001, init, iter.max=10,
	method= c("efron", "breslow", "exact"),
	singular.ok =T, robust=F,
	model=F, x=F, y=T) {

    method <- match.arg(method)
    call <- match.call()
    m <- match.call(expand=F)
    m$method <- m$model <- m$x <- m$y <- m$... <-  NULL
    m$eps <- m$init <- m$iter.max <- m$robust <- m$singular.ok <- NULL
### <TSL>
###    Terms <- if(missing(data)) terms(formula, c('strata', 'cluster'))
###	     else     terms(formula, c('strata', 'cluster'),data=data)    
    Terms <- terms(formula, c('strata', 'cluster'))
### </TSL>
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
###<TSL> we will need this later when dropping terms from formula
    mcopy<-m
###</TSL>
    m <- eval(m, sys.frame(sys.parent()))

    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    weights <- model.extract(m, 'weights')
    offset<- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if(tt == 0)
		    rep(0, nrow(Y))
	      else if(tt == 1)
		      m[[offset]]
	      else {
		    ff <- m[[offset[1]]]
		    for(i in 2:tt)
			    ff <- ff + m[[offset[i]]]
		    ff
		    }

    attr(Terms,"intercept")<- 1  #Cox model always has \Lambda_0
    strats <- attr(Terms, "specials")$strata
    cluster<- attr(Terms, "specials")$cluster
    dropx <- NULL
    if (length(cluster)) {
	if (missing(robust)) robust <- T
	tempc <- untangle.specials(Terms, 'cluster', 1:10)
	ord <- attr(Terms, 'order')[tempc$terms]
	if (any(ord>1)) stop ("Cluster can not be used in an interaction")
	cluster <- strata(m[,tempc$vars], shortlabel=T)  #allow multiples
	dropx <- tempc$terms
	}
    if (length(strats)) {
	temp <- untangle.specials(Terms, 'strata', 1)
	dropx <- c(dropx, temp$terms)
	if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
	else strata.keep <- strata(m[,temp$vars], shortlabel=T)
        ##as.numeric() to codes()
	strats <- codes(strata.keep)
	}
    if (length(dropx)){ 
###<TSL> drop strata/cluster terms from the formula. 
      ##newTerms<-Terms[-dropx]
      ##mcopy$formula<-newTerms
      ##newmlong <- eval(mcopy, sys.frame(sys.parent()))
      ##newm<-m[,names(newmlong),drop=F]
      ##X <- model.matrix(newTerms, newm)[,-1,drop=F]
      X <- model.matrix(Terms[-dropx],m)[,-1,drop=F]
###</TSL>
    }
    else               
      X <- model.matrix(Terms, m)[,-1,drop=F]

    type <- attr(Y, "type")
    if( method=="breslow" || method =="efron") {
	if (type== 'right')  fitter <- get("coxph.fit")
	else if (type=='counting') fitter <- get("agreg.fit")
	else stop(paste("Cox model doesn't support \"", type,
			  "\" survival data", sep=''))
	}
    else if (method=='exact') fitter <- get("agexact.fit")
    else stop(paste ("Unknown method", method))

    if (missing(init)) init <- NULL
    fit <- fitter(X, Y, strats, offset, init=init, iter.max=iter.max,
			eps=eps, weights=weights,
			method=method, row.names(m))

    if (is.character(fit)) {
	fit <- list(fail=fit)
	attr(fit, 'class') <- 'coxph'
	}
    else {
	if (any(is.na(fit$coef))) {
	   vars <- (1:length(fit$coef))[is.na(fit$coef)]
	   msg <-paste("X matrix deemed to be singular; variable",
			   paste(vars, collapse=" "))
	   if (singular.ok) warning(msg)
	   else             stop(msg)
	   }
	fit$n <- nrow(Y)
	attr(fit, "class") <-  fit$method
###<TSL> have to fix up the class of Terms, since it randomly vanishes.
	class(Terms)<-c("terms","formula")
###</TSL>
	fit$terms <- Terms
	fit$assign <- attr(X, 'assign')
	if (robust) {
	    fit$naive.var <- fit$var
	    fit$method    <- method
	    # a little sneaky here: by calling resid before adding the
	    #   na.action method, I avoid having missings re-inserted
	    # I also make sure that it doesn't have to reconstruct X and Y
	    fit2 <- c(fit, list(x=X, y=Y, weights=weights))
	    if (length(strats)) fit2$strata <- strata.keep
	    if (length(cluster)) {
		temp <- residuals.coxph(fit2, type='dfbeta', collapse=cluster,
					  weighted=T)
		indx <- match(cluster, unique(cluster))
		k    <- as.vector(table(indx))
		if (any(k>1)) {
		    #compute the ICC for the m-residuals
		    N    <- sum(k * (k-1))
		    mu   <- sum(fit$resid * (k-1)[indx])/N   #grand mean
		    mu2  <- tapply(fit$resid, indx, mean)    # group means
		    sig  <- tapply((fit$resid - mu)^2, indx, sum)  #SS
		    icc1 <- sum( (k*(mu2-mu))^2 - sig) / sum((k-1)*sig)
		    #rank residuals
		    rr <- rank(fit$resid)
		    mu   <- sum(rr * (k-1)[indx])/N   #grand mean
		    mu2  <- tapply(rr, indx, mean)    # group means
		    sig  <- tapply((rr - mu)^2, indx, sum)  #SS
		    icc2 <- sum( (k*(mu2-mu))^2 - sig) / sum((k-1)*sig)

		    fit$icc <- c(length(k), icc1, icc2)
		    names(fit$icc) <- c("nclust", "icc(resid)",
						    "icc(rank(resid))")
		    }
		}
	    else temp <- residuals.coxph(fit2, type='dfbeta', weighted=T)
	    fit$var <- t(temp) %*% temp
	    }
	na.action <- attr(m, "na.action")
	if (length(na.action)) fit$na.action <- na.action
	if (model) fit$model <- m
	else {
	    if (x)  {
		fit$x <- X
		if (length(strats)) fit$strata <- strata.keep
		}
##<TSL>	    if (y)     fit$y <- Y
	    if (y) {class(Y)<-"Surv";    fit$y <- Y}
##</TSL> (somehow it loses its class)
	    }
	}
    if (!is.null(weights) && any(weights!=1)) fit$weights <- weights

    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- call
    fit$method <- method
    fit
    }

# SCCS @(#)format.Surv.s	4.4 07/10/96
#
# These two functions operate with the newer data.frame code, found on statlib
#   (Eventually part of S, I assume)
#
format.Surv <- function(x, ...) format(as.character.Surv(x), ...)

# The better definition for this is
#   "as.data.frame.Surv <- as.data.frame.model.matrix"
# but, since not everyone has installed the new data.frame software yet, this
# will fail when it can't find as.data.frame.model.matrix.
# So, for this release of survival, the function below is an exact copy of
# as.data.frame.model.matrix, as found on statlib 9/94.
#
#  Changed 5/30/95:  there is a bug in the code I copied: if there are no
#     row names then row.names == character(0), not NULL
as.data.frame.Surv <-
    function(x, row.names = NULL, optional = F)
    {
	d <- dim(x)
	nrows <- d[[1]]
	dn <- dimnames(x)
	row.names <- dn[[1]]
	value <- list(x)
#       if(!is.null(row.names)) {
	if(length(row.names)>0) {
		row.names <- as.character(row.names)
		if(length(row.names) != nrows)
			stop(paste("supplied", length(row.names), 
				"names for a data frame with", nrows, "rows"))
	}
	else if(optional)
		row.names <- character(nrows)
	else row.names <- as.character(seq(length = nrows))
	if(!optional)
		names(value) <- deparse(substitute(x))[[1]]
	attr(value, "row.names") <- row.names
	class(value) <- "data.frame"
	value
}
#
# SCCS @(#)is.ratetable.s	4.5 04/17/96
#
is.ratetable <- function(x, verbose=F) {
    if (!verbose) {
	if (!inherits(x, 'ratetable')) return(F)
	att <- attributes(x)
	if (any(is.na(match(c("dim", "dimnames", "dimid", "factor", "cutpoints"),
			    names(att))))) return(F)
	nd <- length(att$dim)
	if (length(x) != prod(att$dim)) return(F)
	if (!(is.list(att$dimnames) && is.list(att$cutpoints)))
		 return(F)
	if (length(att$dimnames)!=nd || length(att$factor)!=nd ||
			 length(att$cutpoints)!=nd) return(F)
        ## change as.numeric() to codes()
	fac <- as.numeric(att$factor)
	if (any(is.na(fac))) return(F)
	if (any(fac <0)) return(F)
	for (i in 1:nd) {
	    n <- att$dim[i]
	    if (length(att$dimnames[[i]]) !=n) return(F)
	    if (fac[i]!=1 && length(att$cutpoints[[i]])!=n) return(F)
	    if (fac[i]!=1 && any(order(att$cutpoints[[i]])!= 1:n)) return(F)
	    if (fac[i]==1 && !is.null(att$cutpoints[[i]]))  return(F)
	    if (fac[i]>1 && i<nd) return(F)
	    }
	return(T)
	}

    #verbose return messages, useful for debugging
    msg <- ""
    if (!inherits(x, 'ratetable')) msg <- c(msg, "wrong class")
    att <- attributes(x)
    if (any(is.na(match(c("dim", "dimnames", "dimid", "factor", "cutpoints"),
			names(att))))) msg <- c(msg, 'missing an attribute')
    nd <- length(att$dim)
    if (length(x) != prod(att$dim)) msg <- c(msg, 'dims dont match length')
    if (!(is.list(att$dimnames) && is.list(att$cutpoints)))
	     msg <- c(msg, 'dimnames or cutpoints not a list')
    if (length(att$dimnames)!=nd || length(att$factor)!=nd ||
		     length(att$cutpoints)!=nd) msg <- c(msg, 'bad lengths')
    ##change as.numeric() to codes()
    fac <- as.numeric(att$factor)
    if (any(is.na(fac))) msg <- c(msg, 'a missing factor')
    if (any(fac <0)) msg <- c(msg, 'factor <0')
    for (i in 1:nd) {
	n <- att$dim[i]
	if (length(att$dimnames[[i]]) !=n) msg <- c(msg, 'dimnames wrong length')
	if (fac[i]!=1 && length(att$cutpoints[[i]])!=n) msg <- c(msg, 'cutpnt missing')
	if (fac[i]!=1 && any(order(att$cutpoints[[i]])!= 1:n)) msg <- c(msg, 'unsorted cutpoints')
	if (fac[i]==1 && !is.null(att$cutpoints[[i]]))  msg <- c(msg, 'cutpnt should be null')
	if (fac[i]>1 && i<nd) msg <- c(msg, 'only the last dim can be interpolated')
	}
    if (msg=='') T
    else msg
    }
#
# SCCS @(#)is.ratetable.verbose.s	1.1 10/27/95
#   A version of the function that tells you WHY it's not a ratetable
#

is.ratetable.verbose <- function(x) {
    if (!inherits(x, 'ratetable')) return("wrong class")
    att <- attributes(x)
    if (any(is.na(match(c("dim", "dimnames", "dimid", "factor", "cutpoints"),
			names(att))))) return('missing an attribute')
    nd <- length(att$dim)
    if (length(x) != prod(att$dim)) return('dims dont match length')
    if (!(is.list(att$dimnames) && is.list(att$cutpoints)))
	     return('dimnames or cutpoints not a list')
    if (length(att$dimnames)!=nd || length(att$factor)!=nd ||
		     length(att$cutpoints)!=nd) return('bad lengths')
    ##as.numeric() to codes()
    fac <- as.numeric(att$factor)
    if (any(is.na(fac))) return('a missing factor')
    if (any(fac <0)) return('factor <0')
    for (i in 1:nd) {
	n <- att$dim[i]
	if (length(att$dimnames[[i]]) !=n) return('dimnames wrong length')
	if (fac[i]!=1 && length(att$cutpoints[[i]])!=n) return('cutpnt missing')
	if (fac[i]!=1 && any(order(att$cutpoints[[i]])!= 1:n)) return('unsorted cutpoints')
	if (fac[i]==1 && !is.null(att$cutpoints[[i]]))  return('cutpnt should be null')
	if (fac[i]>1 && i<nd) return('only the last dim can be interpolated')
	}
    T
    }
# SCCS @(#)lines.survfit.s	4.8  12/14/94
lines.survfit <- function(x, type='s', mark=3, col=1, lty=1, lwd=1,
		       mark.time =T, xscale=1, yscale=1,  ...) {
    if (inherits(x, 'survexp')) {
	if (missing(type)) type <- 'l'
	if (!is.numeric(mark.time)) mark.time <- F
	}
    if (is.numeric(mark.time)) mark.time <- sort(unique(mark.time[mark.time>0]))

    if (is.null(x$strata)) {
	if (is.matrix(x$surv)) ncurv <- ncol(x$surv)
	else ncurve <- 1
	nstrat <- 1
	strata <- length(x$surv)/nstrat
	}
    else {
	strata <- x$strata
	nstrat <- length(strata)
	ncurve <- nstrat * length(x$surv)/ sum(strata)
	}

    mark <- rep(mark, length=ncurve)
    col  <- rep(col , length=ncurve)
    lty  <- rep(lty , length=ncurve)
    lwd  <- rep(lwd , length=ncurve)
    time <- rep(x$time, length=length(x$surv))
    j <- 1
    for (i in 1:ncurve) {
	n <- strata[1+(i-1)%%nstrat]
	who <- seq(from=j, length=n)
	j <-  j+n
	xx <- c(0, time[who])/xscale
	yy <- c(1, x$surv[who])*yscale
	nn <- length(xx)
	#
	# This 'whom' addition is to replace verbose horizonal sequences
	#  like (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
	#  with (1, .2), (3, .1) if type='s' and (1, .2), (2.9, .2), (3, .1)
	#  otherwise.  They are slow, and can smear the looks of a line type
	#
	whom <- c(match(unique(yy[-nn]), yy), nn)
	if (type!='s') whom <- sort(unique(c(whom, whom[-1]-1)))
	lines(xx[whom], yy[whom], type=type, col=col[i], lty=lty[i], lwd=lwd[i], ...)

	if (is.numeric(mark.time)) {
	    indx <- mark.time
	    for (k in seq(along=mark.time))
		indx[k] <- sum(mark.time[k] > xx)
	    points(mark.time[indx<nn], yy[indx[indx<nn]],
		   pch=mark[i],col=col[i], ...)
	    }
	else if (mark.time==T) {
	    deaths <- c(-1, x$n.event[who])
	    if ( any(deaths==0))
		points(xx[deaths==0], yy[deaths==0],
			      pch=mark[i],col=col[i], ...)
	    }
	}
    invisible()
    }
# SCCS @(#)match.ratetable.s	4.4 05/31/94
# Do a set of error checks on whether the ratetable() vars match the
#   actual ratetable
# This is called by pyears and survexp, but not by users
#
# Returns a subscripting vector and a call
#
match.ratetable <- function(R, ratetable) {
    attR <- attributes(R)
    attributes(R) <- attR['dim']     #other attrs get in the way later
    if (!is.ratetable(ratetable)) stop("Invalid rate table")
    dimid <- attr(ratetable, 'dimid')
    ord <- match(attR$dimnames[[2]], dimid)
    if (any(is.na(ord)))
       stop(paste("Argument '", (attR$dimnames[[2]])[is.na(ord)],
	    "' in ratetable()",
	    " does not match the given table of event rates", sep=''))
    nd <- length(ord)
    if (nd != length(dimid))
	stop("The ratetable() call has the wrong number of arguments")
    ##<TSL> ord[ord]<- doesn't work, but I can't reproduce the bug
    ord1<-ord
    ord[ord1] <- 1:nd   #reverse the index, so "ord" can be on the right-hand
    ##</TSL>
    R <- R[,ord,drop=F]

    # Check out the dimensions of R --
    const <- attR[["constants"]][ord]
    call <- "ratetable["
    levlist <- attR[['levlist']][ord]
    dtemp <-dimnames(ratetable)
    efac  <- attr(ratetable, 'factor')
    for (i in (1:nd)) {
	if (const[i]) {   #user put in a constant
	    temp <- match(levlist[[i]], dtemp[[i]])
	    if (is.na(temp)) {
                ## as.numeric() to codes()
		temp <- as.numeric(levlist[[i]])
		if (is.na(temp))
		       stop(paste("Invalid value in ratetable() for variable",
				 dimid[i]))
		if (efac[i]==1) {  # this level is a factor
		    if (temp<=0 || temp!=floor(temp) || temp >length(dtemp[[i]]))
		       stop(paste("Invalid value in ratetable() for variable",
				 dimid[i]))
		    }
		else stop(paste("Invalid value in ratetable() for variable",
					dimid[i]))
		}
	    R[,i] <- temp
	    call <- paste(call, temp)
	    }
	else if (length(levlist[[i]]) >0) {  #factor or character variable
	    if (efac[i]!=1) stop(paste("In ratetable(),", dimid[i],
				     "must be a continuous variable"))
	    temp <- match(levlist[[i]], dtemp[[i]])
	    if (any(is.na(temp)))
		stop(paste("Levels do not match for ratetable() variable",
			    dimid[i]))
	    R[,i] <- temp[R[,i]]
	    }
	else {   # ratetable() thinks it is a continuous variable
	    if (efac[i]==1) {   #but it's not-- make sure it is an integer
		temp <- R[,i]
		if (any(floor(temp)!=temp) || any(temp<0) ||
			    max(temp) > length(dtemp[[i]]))
		stop(paste("In ratetable()",dimid[i],"out of range"))
		}
	    }
	if (i==nd) call <- paste(call, "]")
	else       call <- paste(call, ",")
	}

    summ <- attr(ratetable, 'summary')
    if (is.null(summ))
	 list(R= R[,!const, drop=F], call={if(any(const)) call else NULL})
    else list(R= R[,!const, drop=F], call={if(any(const)) call else NULL},
		summ=summ(R))
    }
model.frame.coxph <- function(object, ...) {
    Call <- object$call
    Call[[1]] <- as.name("model.frame")
    Call <- Call[match(c("", "formula", "data", "weights", "subset",
			   "na.action"), names(Call), 0)]
    eval(as.call(Call))
    }
#SCCS 06/10/92 @(#)model.frame.default.s	4.3
# Only change -- look to options() for the default na.action
#


naprint.omit <- function(x)
    paste(length(x), "deleted due to missing")

# Put the missing values back into a vector.
#   And be careful about the labels too.
naresid.omit <- function(omit, x) {
    if (length(omit)==0 || !is.numeric(omit))
	stop("Invalid argument for 'omit'")

    if (is.matrix(x)) {
	n <- nrow(x)
	keep <- rep(NA,n+ length(omit))
	keep[-omit] <- 1:n
	x <- x[keep,,drop=F]
	temp <- dimnames(x)[[1]]
	if (length(temp)) {
	    temp[omit] <- names(omit)
	    dimnames(x) <- list(temp, dimnames(x)[[2]])
	    }
	x
	}
    else {
	n <- length(x)
	keep <- rep(NA,n+ length(omit))
	keep[-omit] <- 1:n
	x <- x[keep]
	temp <- names(x)
	if (length(temp)) {
	    temp[omit] <- names(omit)
	    names(x) <- temp
	    }
	x
	}
    }
naprint.omit <- function(x)
    paste(length(x), "observations deleted due to missing")
naprint <- function(x, ...)
    UseMethod("naprint",x,...)

naprint.default <- function(...)
    return("")
#SCCS 04/14/92 @(#)naresid.omit.s	4.2
naresid.omit <- function(omit, x) {
    if (length(omit)==0 || !is.numeric(omit))
	stop("Invalid argument for 'omit'")

    if (is.matrix(x)) {
	n <- nrow(x)
	keep <- rep(NA,n+ length(omit))
	keep[-omit] <- 1:n
	x <- x[keep,,drop=F]
	temp <- dimnames(x)[[1]]
	if (length(temp)) {
	    temp[omit] <- names(omit)
	    dimnames(x) <- list(temp, dimnames(x)[[2]])
	    }
	x
	}
    else {
	n <- length(x)
	keep <- rep(NA,n+ length(omit))
	keep[-omit] <- 1:n
	x <- x[keep]
	temp <- names(x)
	if (length(temp)) {
	    temp[omit] <- names(omit)
	    names(x) <- temp
	    }
	x
	}
    }
naresid <- function(omit, ...)
	UseMethod("naresid",omit,...)

naresid.default <- function(omit, x)  x
#SCCS @(#)plot.cox.zph.s	4.6 08/13/96
plot.cox.zph <- function(x, resid=T, se=T, df=4, nsmo=40, var, ...) {
    require("splines")
    xx <- x$x
    yy <- x$y
    d <- nrow(yy)
    df <- max(df)     #error proofing
    nvar <- ncol(yy)
    pred.x <- seq(from=min(xx), to=max(xx), length=nsmo)
    temp <- c(pred.x, xx)
    lmat <- ns(temp, df=df, intercept=T)
    pmat <- lmat[1:nsmo,]       # for prediction
    xmat <- lmat[-(1:nsmo),]
    qmat <- qr(xmat)

    if (se) {
	bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
	xtx <- bk %*% t(bk)
	seval <- d*((pmat%*% xtx) *pmat) %*% rep(1, df)
	}

    ylab <- paste("Beta(t) for", dimnames(yy)[[2]])
    if (missing(var)) var <- 1:nvar
    else {
	if (is.character(var)) var <- match(var, dimnames(yy)[[2]])
	if  (any(is.na(var)) || max(var)>nvar || min(var) <1)
	    stop("Invalid variable requested")
	}

    #
    # Figure out a 'good' set of x-axis labels.  Find 8 equally spaced
    #    values on the 'transformed' axis.  Then adjust until they correspond
    #    to rounded 'true time' values.  Avoid the edges of the x axis, or
    #    approx() may give a missing value
    if (x$transform == 'log') {
	xx <- exp(xx)
	pred.x <- exp(pred.x)
	}
    else if (x$transform != 'identity') {
	xtime <- as.numeric(dimnames(yy)[[1]])
	apr1  <- approx(xx, xtime, seq(min(xx), max(xx), length=17)[2*(1:8)])
	temp <- signif(apr1$y,2)
	apr2  <- approx(xtime, xx, temp)
	xaxisval <- apr2$y
	xaxislab <- rep("",8)
	for (i in 1:8) xaxislab[i] <- format(temp[i])
	}

    for (i in var) {
	y <- yy[,i]
	yhat <- pmat %*% qr.coef(qmat, y)
	if (resid) yr <-range(yhat, y)
	else       yr <-range(yhat)
	if (se) {
	    temp <- 2* sqrt(x$var[i,i]*seval)
	    yup <- yhat + temp
	    ylow<- yhat - temp
	    yr <- range(yr, yup, ylow)
	    }

	if (x$transform=='identity')
	    plot(range(xx), yr, type='n', xlab="Time", ylab=ylab[i], ...)
	else if (x$transform=='log')
	    plot(range(xx), yr, type='n', xlab="Time", ylab=ylab[i], log='x',
			...)
	else {
	    plot(range(xx), yr, type='n', xlab="Time", ylab=ylab[i], axes=F,...)
	    axis(1, xaxisval, xaxislab)
	    axis(2)
	    box()
	    }
	if (resid) points(xx, y)
	lines(pred.x, yhat)
	if (se) {
	    lines(pred.x, yup,lty=2)
	    lines(pred.x, ylow, lty=2)
	    }
	}
    }
#SCCS 04/14/92 @(#)plot.coxph.s	4.2
plot.coxph <- function(fit, ...) {
    op <- par(ask=T)
    on.exit(par(op))
    yy <- (1-fit$residuals)/ fit$linear.predictors   # psuedo y
    plot(fit$linear.predictors, rank(yy))

    std.resid <- fit$residuals/ sqrt(predict(fit, type='expected'))
    plot(fit$linear.predictors, std.resid,
	xlab='Linear predictor', ylab='standardized residual')

    }
#SCCS @(#)plot.survfit.s	4.11 12/14/94
plot.survfit<- function(surv, conf.int,  mark.time=T,
		 mark=3,col=1,lty=1, lwd=1, cex=1,log=F, yscale=1,
		 xscale=1,
		 xlab="", ylab="", xaxs='i', ...) {

    if (!inherits(surv, 'survfit'))
	  stop("First arg must be the result of survfit")

    stime <- surv$time / xscale
    ssurv <- surv$surv
    if (missing(conf.int)) {
	if (is.null(surv$strata) && !is.matrix(ssurv)) conf.int <-T
	else conf.int <- F
	}

    if (is.null(surv$strata)) {
	nstrat <- 1
	stemp <- rep(1, length(surv$time))
	}
    else {
	nstrat <- length(surv$strata)
	stemp <- rep(1:nstrat,surv$strata)
	}
    if (is.null(surv$n.event)) mark.time <- F   #expected survival curve

    # set default values for missing parameters
    if (is.matrix(ssurv)) ncurve <- nstrat * ncol(ssurv)
    else                  ncurve <- nstrat
    mark <- rep(mark, length=ncurve)
    col  <- rep(col, length=ncurve)
    lty  <- rep(lty, length=ncurve)
    lwd  <- rep(lwd, length=ncurve)
###<TSL>
###    if (is.numeric(mark.time)) mark.time <- sort(mark.time[mark.time>0])
    if (!is.logical(mark.time) && is.numeric(mark.time)) mark.time <- sort(mark.time[mark.time>0])
###</TSL>
    if (missing(xaxs)) temp <- 1.04*max(stime)
    else               temp <- max(stime)
    #
    # for log plots we have to be tricky about the y axis scaling
    #
    if   (log) {
	    ymin <- min(.1,ssurv[!is.na(ssurv) &ssurv>0])
	    ssurv[!is.na(ssurv) &ssurv==0] <- ymin
	    plot(c(0, temp),
	       yscale*c(.99, ymin),
	       type ='n', log='y', xlab=xlab, ylab=ylab, xaxs=xaxs,...)
	    }
     else
	 plot(c(0, temp), yscale*c(0,1),
	      type='n', xlab=xlab, ylab=ylab, xaxs=xaxs, ...)

    if (yscale !=1) par(usr=par("usr")/ c(1,1,yscale, yscale))
    #
    # put up the curves one by one
    #   survfit has already put them into the "right" order
    i _ 0
    xend _ NULL
    yend _ NULL

    #
    # The 'whom' addition is to replace verbose horizonal sequences
    #  like (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
    #  with (1, .2), (3, .1) -- remember that type='s'.  They are slow,
    #  and can smear the looks of a line type
    #
    for (j in unique(stemp)) {
	who _ (stemp==j)
	xx _ c(0,stime[who])
	deaths <- c(-1, surv$n.event[who])
	if (is.matrix(ssurv)) {
	    for (k in 1:ncol(ssurv)) {
		i _ i+1
		yy _ c(1,ssurv[who,k])
		nn <- length(xx)
		whom <- c(match(unique(yy[-nn]), yy), nn)
		lines(xx[whom], yy[whom], lty=lty[i], col=col[i], lwd=lwd[i], type='s')

		if (is.numeric(mark.time)) {
		    indx <- mark.time
		    for (k in seq(along=mark.time))
			indx[k] <- sum(mark.time[k] > xx)
		    points(mark.time[indx<nn], yy[indx[indx<nn]],
			   pch=mark[i],col=col[i],cex=cex)
		    }
		else if (mark.time==T && any(deaths==0))
		    points(xx[deaths==0], yy[deaths==0],
			   pch=mark[i],col=col[i],cex=cex)
		xend _ c(xend,max(xx))
		yend _ c(yend,min(yy))

		if (conf.int==T && !is.null(surv$upper)) {
		    if (ncurve==1) lty[i] <- lty[i] +1
		    yy _ c(1,surv$upper[who,k])
		    lines(xx,yy, lty=lty[i], col=col[i], lwd=lwd[i], type='s')
		    yy _ c(1,surv$lower[who,k])
		    lines(xx,yy, lty=lty[i], col=col[i], lwd=lwd[i], type='s')
		    }
		}
	    }

	else {
	    i <- i+1
	    yy _ c(1,ssurv[who])
	    nn <- length(xx)
	    whom <- c(match(unique(yy[-nn]), yy), nn)
	    lines(xx[whom], yy[whom], lty=lty[i], col=col[i], lwd=lwd[i], type='s')

	    if (!is.logical(mark.time) && is.numeric(mark.time)) {
		indx <- mark.time
		for (k in seq(along=mark.time))
		    indx[k] <- sum(mark.time[k] > xx)
		points(mark.time[indx<nn], yy[indx[indx<nn]],
		       pch=mark[i],col=col[i],cex=cex)
		}
	    else if (mark.time==T && any(deaths==0))
		points(xx[deaths==0], yy[deaths==0],
		       pch=mark[i],col=col[i],cex=cex)
	    xend _ c(xend,max(xx))
	    yend _ c(yend,min(yy))

	    if (conf.int==T && !is.null(surv$upper)) {
		if (ncurve==1) lty[i] <- lty[i] +1
		yy _ c(1,surv$upper[who])
		lines(xx,yy, lty=lty[i], col=col[i], lwd=lwd[i], type='s')
		yy _ c(1,surv$lower[who])
		lines(xx,yy, lty=lty[i], col=col[i], lwd=lwd[i], type='s')
		}
	    }
	}
    invisible(list(x=xend, y=yend))
}

#SCCS 04/14/92 @(#)points.survfit.s	4.2
points.survfit <- function(object, ...) {
    if (!is.matrix(object$surv))
	    points(object$time, object$surv, ...)
    else
	    matpoints(object$time, object$surv, ...)
    }
#SCCS 06/29/93 @(#)predict.coxph.s	4.8
#What do I need to do predictions --
#
#linear predictor:  exists
#        +se     :  X matrix
#        +newdata:  means of old X matrix, new X matrix, new offset
#
#risk -- same as lp
#
#expected --    cumulative hazard for subject= baseline haz + time + risk
#        +se :  sqrt(expected)
#      +new  :  baseline hazard function, new time, new x, means of old X,
#                        new offset, new strata
#
#terms -- : X matrix and the means
#    +se  :  ""  + I matrix
#   +new  : new X matrix and the old means + I matrix
predict.coxph <-
function(object, newdata, type=c("lp", "risk", "expected", "terms"),
		se.fit=F,
		terms=labels.lm(object), collapse, safe=F, ...)

    {
    type <-match.arg(type)
    n <- object$n
    Terms <- object$terms
    strata <- attr(Terms, 'specials')$strata
    if (length(strata)) {
	   temp <- untangle.specials(Terms, 'strata', 1)
	   Terms2 <- Terms[-temp$terms]
	   }
    else  Terms2 <- Terms
    offset <- attr(Terms, "offset")
### <TSL> this is never used, so it doesn't matter that it doesn't work
###    resp <- attr(Terms, "variables")[attr(Terms, "response")]

    if (missing(newdata)) {
	if (type=='terms' || (se.fit && (type=='lp' || type=='risk'))) {
	    x <- object$x
	    if (is.null(x)) {
		x <- model.matrix(Terms2, model.frame(object))[,-1,drop=F]
		}
###	    x <- sweep(x, 2, object$means)
	    x<-x-object$means
	    }
	else if (type=='expected') {
	    y <- object$y
	    if (missing(y)) {
		m <- model.frame(object)
		y <- model.extract(m, 'response')
		}
	    }
	}
    else {
	if (type=='expected')
	     m <- model.newframe(Terms, newdata, response=T)
	else m <- model.newframe(Terms2, newdata)

	x <- model.matrix(Terms2, m)[,-1,drop=F]
###	x <- sweep(x, 2, object$means)
	x<- x -object$means
	if (length(offset)) {
	    if (type=='expected') offset <- as.numeric(m[[offset]])
	    else {
		offset <- attr(Terms2, 'offset')
		offset <- as.numeric(m[[offset]])
		}
	    }
	else offset <- 0
	}

    #
    # Now, lay out the code one case at a time.
    #  There is some repetition this way, but otherwise the code just gets
    #    too complicated.
    coef <- ifelse(is.na(object$coef), 0, object$coef)
    if (type=='lp' || type=='risk') {
	if (missing(newdata)) {
	    pred <- object$linear.predictors
	    names(pred) <- names(object$residuals)
	    }
	else                  pred <- x %*% coef  + offset
	if (se.fit) se <- sqrt(diag(x %*% object$var %*% t(x)))

	if (type=='risk') {
	    pred <- exp(pred)
	    if (se.fit) se <- se * sqrt(pred)
	    }
	}

    else if (type=='expected') {
	if (missing(newdata)) pred <- y[,ncol(y)] - object$residual
	else  stop("Method not yet finished")
	se   <- sqrt(pred)
	}

    else {  #terms
	attr(x, "constant") <- rep(0, ncol(x))
	asgn <- object$assign
	terms <- match.arg(Terms2, labels.lm(object))
	asgn <- asgn[terms]
	if (se.fit) {
	    temp <- Build.terms(x, coef, object$var, asgn, F)
	    pred <- temp$fit
	    se   <- temp$se.fit
	    }
	else pred<- Build.terms(x, coef, NULL, asgn, F)
	}

    if (se.fit) se <- drop(se)
    pred <- drop(pred)
    #Expand out the missing values in the result
    if (!is.null(object$na.action)) {
	pred <- naresid(object$na.action, pred)
	if (is.matrix(pred)) n <- nrow(pred)
	else               n <- length(pred)
	if(se.fit) se <- naresid(object$na.action, se)
	}

    # Collapse over subjects, if requested
    if (!missing(collapse)) {
	if (length(collapse) != n) stop("Collapse vector is the wrong length")
	pred <- rowsum(pred, collapse)
	if (se.fit) se <- sqrt(rowsum(se^2, collapse))
	}

    if (se.fit) list(fit=pred, se.fit=se)
    else pred
    }
# SCCS @(#)predict.survreg.s	4.7 07/14/92
predict.survreg <-
    function(object, newdata, ...)
    {
    # This routine really hasn't been written-- you are looking at a
    #          placeholder.  Some glm aspects do work, however.
    NextMethod("predict")
    }
# SCCS @(#)print.cox.zph.s	4.5 09/27/96
print.cox.zph <- function(x, digits = max(options()$digits - 4, 3))
    invisible(print(x$table, digits=digits))
#SCCS 09/27/96 @(#)print.coxph.null.s	4.6
print.coxph.null <-
 function(cox, digits=max(options()$digits - 4, 3), ...)
    {
    if (!is.null(cl<- cox$call)) {
	cat("Call:  ")
	dput(cl)
	cat("\n")
	}

    cat("Null model\n  log likelihood=", format(cox$loglik), "\n")
    omit <- cox$na.action
    if (length(omit))
	cat("  n=", cox$n, " (", naprint(omit), ")\n",
				sep="")
    else cat("  n=", cox$n, "\n")
    }
# SCCS @(#)print.coxph.s	4.8 07/10/96
print.coxph <-
 function(cox, digits=max(options()$digits - 4, 3), ...)
    {
    if (!is.null(cl<- cox$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
    if (!is.null(cox$fail)) {
	cat(" Coxph failed.", cox$fail, "\n")
	return()
	}
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    coef <- cox$coef
    se <- sqrt(diag(cox$var))
    if(is.null(coef) | is.null(se))
        stop("Input is not valid")

    if (is.null(cox$naive.var)) {
	tmp <- cbind(coef, exp(coef), se, coef/se,
	       signif(1 - pchisq((coef/ se)^2, 1), digits -1))
	dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
	    "se(coef)", "z", "p"))
	}
    else {
	nse <- sqrt(diag(cox$naive.var))
	tmp <- cbind(coef, exp(coef), nse, se, coef/se,
	       signif(1 - pchisq((coef/se)^2, 1), digits -1))
	dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
	    "se(coef)", "robust se", "z", "p"))
	}
    cat("\n")
    surv4.prmatrix(tmp)

    logtest <- -2 * (cox$loglik[1] - cox$loglik[2])
    df <- sum(!is.na(coef))
    cat("\n")
    cat("Likelihood ratio test=", format(round(logtest, 2)), "  on ",
	df, " df,", " p=", format(1 - pchisq(logtest, df)),  sep="")
    omit <- cox$na.action
    if (length(omit))
	cat("  n=", cox$n, " (", naprint(omit), ")\n", sep="")
    else cat("  n=", cox$n, "\n")
    if (length(cox$icc))
	cat("   number of clusters=", cox$icc[1],
	    "    ICC=", format(cox$icc[2:3]), "\n")
    invisible()
    }
#SCCS 01/03/94 @(#)print.data.frame.s	4.4
#
# In order to get objects with attributes to print correctly, I replace the
#   call to "as.matrix" with a copy of as.matrix.data.frame, one that knows
#   its output is character, and so calls the appropriate as.character routine
#
##<TSL> disable this by renaming - it doesn't work </TSL>
dont.print.data.frame <-
function(x, ..., quote = F, right = T)
{
    if(length(x)==0)
        cat("NULL data frame with", length(row.names(x)), "rows\n")

    else if(length(row.names(x))==0) {
        print.atomic(names(x), quote = F)
        cat("< 0 rows> \n")
        }
    else {
	DD <- dim(x)
	dn <- dimnames(x)
	collabs <- as.list(dn[[2]])
	class(x) <- NULL
	p <- DD[2]
	n <- DD[1]
	non.numeric <- non.atomic <- F
	for(j in 1:p) {
	    xj <- x[[j]]

	    # 4 Line below are the ones I added
	    if (!is.null(cl <-class(xj)) &&
		       (cl == 'Surv' || cl=='date'))
			     xj <- x[[j]] <- as.character(x[[j]])
	    if(length(dj <- dim(xj)) == 2 && dj[2] > 1) {
		if (inherits(xj, "data.frame"))
		    xj <- x[[j]] <- as.matrix(x[[j]])
		dnj <- dimnames(xj)[[2]]
		collabs[[j]] <- paste(collabs[[j]], if(length(dnj)) dnj
		     else seq(1:dj[2]), sep = ".")
		}
	    if(length(levels(xj)) > 0 || (!is.complex(xj) ||!is.numeric(xj)))
		non.numeric <- TRUE
	    if(!is.atomic(xj))
		non.atomic <- TRUE
	    }
	if(non.atomic) {
	    for(j in 1:p) {
		xj <- x[[j]]
		if(is.recursive(xj)) {
		    }
		else x[[j]] <- as.list(as.vector(xj))
		}
	    }
	else if(non.numeric) {
	    for(j in 1:p) {
		xj <- x[[j]]
		if(length(levels(xj)))
		    x[[j]] <- as.vector(xj)
		else x[[j]] <- format(xj)
		}
	    }
	x <- unlist(x, recursive = F)
	dim(x) <- c(n, length(x)/n)
	dimnames(x) <- list(dn[[1]], unlist(collabs, use.names=F))
	class(x) <- "matrix"
	print(x, ..., quote=quote, right=right)
	}
    invisible(x)
    }
#SCCS @(#)print.summary.survfit.s	4.2 09/27/96
print.summary.survfit <- function(fit, digits = max(options()$digits - 4, 3), ...) {
    savedig <- options(digits=digits)
    on.exit(options(savedig))

    if (!is.null(cl<- fit$call)) {
	cat("Call: ")
	dput(cl)
	cat("\n")
	}

    omit <- fit$na.action
    if (length(omit))
	cat(naprint(omit), "\n")

    mat <- cbind(fit$time, fit$n.risk, fit$n.event, fit$surv)
    cnames <- c("time", "n.risk", "n.event")
    if (is.matrix(fit$surv)) ncurve <- ncol(fit$surv)
    else ncurve <- 1

    if (ncurve==1) {                 #only 1 curve
	cnames <- c(cnames, "survival")
	if (!is.null(fit$std.err)) {
	    if (is.null(fit$lower)) {
		mat <- cbind(mat, fit$std.err)
		cnames <- c(cnames, "std.err")
		}
	    else {
		mat <- cbind(mat, fit$std.err, fit$lower, fit$upper)
		cnames <- c(cnames, 'std.err',
			  paste("lower ", fit$conf.int*100, "% CI", sep=''),
			  paste("upper ", fit$conf.int*100, "% CI", sep=''))
		}
	    }
	}
    else cnames <- c(cnames, paste("survival", seq(ncurve), sep=''))

    dimnames(mat) <- list(NULL, cnames)
    if (is.null(fit$strata)) {
	surv4.prmatrix(mat, rowlab=rep("", nrow(mat)))
	}
    else  { #print it out one strata at a time
	for (i in levels(fit$strata)) {
	    who <- (fit$strata==i)
	    cat("               ", i, "\n")
	    if (sum(who) ==1) print(mat[who,])
	    else    surv4.prmatrix(mat[who,], rowlab=rep("", sum(who)))
	    cat("\n")
	    }
	}
    invisible(fit)
    }
# SCCS @(#)print.summary.survreg.s	4.8 09/27/96
print.summary.survreg <- function(x, digits = max(options()$digits - 4, 3), quote = T, prefix = "")
{
    nas <- x$nas
    coef <- x$coef
    correl <- x$correl
    if(any(nas)) {
        nc <- length(nas)
        cnames <- names(nas)
        coef1 <- array(NA, c(nc, 3), list(cnames, dimnames(coef)[[2]]))
            
        coef1[!nas,  ] <- coef
        coef <- coef1
        if(!is.null(correl)) {
            correl1 <- matrix(NA, nc, nc, dimnames = list(cnames,
                cnames))
            correl1[!nas, !nas] <- correl
            correl <- correl1
            }
        }
    if(is.null(digits))
        digits <- options()$digits
    else options(digits = digits)
    cat("\nCall:\n")
    dput(x$call)
    dresid <- x$deviance.resid
    n <- length(dresid)
    rdf <- x$df[2]
    if(rdf > 5) {
        cat("Deviance Residuals:\n")
        rq <- quantile(as.vector(dresid))
        names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
        print(rq, digits = digits)
        }
    else if(rdf > 0) {
        cat("Deviance Residuals:\n")
        print(dresid, digits = digits)
        }
    if(any(nas))
        cat("\nCoefficients: (", sum(nas), 
            " not defined because of singularities)\n", sep = "")
        
    else cat("\nCoefficients:\n")
    print(coef, digits = digits)
    omit <- x$na.action
    if (length(omit))
	cat("  n=", n, " (", naprint(omit), ")\n", sep="")

    cat("\n", x$parms, "\n", sep='')
    int <- attr(x$terms, "intercept")
    if(is.null(int))
        int <- 1
    temp <- format(round(c(x$null.deviance, x$deviance), digits))
    cat("\n    Null Deviance:", temp[1], "on",
		     n - int, "degrees of freedom\n")
    cat("Residual Deviance:", temp[2], "on",
	   round(rdf, digits), "degrees of freedom  (LL=",
		format(x$loglik), ")\n")
    cat("Number of Newton-Raphson Iterations:", format(trunc(x$iter)),
        "\n")
    if(!is.null(correl)) {
        p <- dim(correl)[2]
        if(p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            ll <- lower.tri(correl)
            correl[ll] <- format(round(correl[ll], digits))
            correl[!ll] <- ""
            print(correl[-1,  - p, drop = F], quote = F, digits = 
                digits)
            }
        }
    cat("\n")
    invisible(NULL)
    }
#SCCS 09/27/96 @(#)print.survdiff.s	4.10
print.survdiff <- function(fit, digits = max(options()$digits - 4, 3), ...) {

    saveopt <-options(digits=digits)
    on.exit(options(saveopt))

    if (!inherits(fit, 'survdiff'))
	stop("Object is not the result of survdiff")
    if (!is.null(cl<- fit$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}

    omit <- fit$na.action
    if (length(omit)) cat("n=", sum(fit$n), ", ", naprint(omit),
					  ".\n\n", sep='')

    if (length(fit$n)==1)  {
	z <- sign(fit$exp - fit$obs) * sqrt(fit$chisq)
	temp <- c(fit$obs, fit$exp, z, signif(1-pchisq(fit$chisq, 1),digits))
	names(temp) <- c("Observed", "Expected", "Z", "p")
	print(temp)
	}
    else {
	if (is.matrix(fit$obs)){
	    otmp <- apply(fit$obs,1,sum)
	    etmp <- apply(fit$exp,1,sum)
	    }
	else {
	    otmp <- fit$obs
	    etmp <- fit$exp
	    }
	df <- (sum(1*(etmp>0))) -1
##	temp <- cbind(fit$n, otmp, etmp, ((otmp-etmp)^2)/ etmp,((otmp-etmp)^2)/ diag(fit$var))
##	dimnames(temp) <- list(names(fit$n), c("N", "Observed", "Expected","(O-E)^2/E", "(O-E)^2/V"))
	temp <- cbind(as.vector(fit$n), otmp, etmp, ((otmp - etmp)^2)/etmp, ((otmp - etmp)^2)/diag(fit$var))
	dimnames(temp) <- list(rownames(fit$n), c("N", "Observed", "Expected", "(O-E)^2/E", "(O-E)^2/V"))
	print(temp)
	cat("\n Chisq=", format(round(fit$chisq,1)),
		 " on", df, "degrees of freedom, p=",
		 format(signif(1-pchisq(fit$chisq, df),digits)), "\n")
       }
    invisible(fit)
    }
#SCCS @(#)print.survexp.s	4.11 09/27/96
print.survexp <- function(fit, scale=1, digits = max(options()$digits - 4, 3), naprint=F, ...) {
    if (!inherits(fit, 'survexp'))
	    stop("Invalid data")
    savedig <- options(digits=digits)
    on.exit(options(savedig))

    if (!is.null(cl<- fit$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}

    if (!is.null(fit$summ)) cat(fit$summ)
    omit <- fit$na.action
    if (length(omit))
	cat(naprint(omit), "\n")
    else cat("\n")

    if (is.null(fit$strata))  { #print it as a matrix
	mat <- cbind(fit$time/scale, fit$n.risk, fit$surv, fit$std.err)
	if (!naprint) {
	    miss <- (is.na(mat)) %*% rep(1,ncol(mat))
	    mat <- mat[miss<(ncol(mat)-2),,drop=F]
	    }
	if (is.matrix(fit$surv)) cname <- dimnames(fit$surv)[[2]]
	else                     cname <- "survival"
	if (!is.null(fit$std.err))
	      cname <- c(cname, paste("se(", cname, ")", sep=''))
	surv4.prmatrix(mat, rowlab=rep("", nrow(mat)),
		   collab=c("Time", "n.risk", cname))
	}
    else  { #print it out one strata at a time, since n's differ
	if (is.null(fit$std.err)) tname <- 'survival'
	else                      tname <- c('survival', 'se(surv)')
	nstrat <- length(fit$strata)
	levs <- names(fit$strata)
	if (nrow(fit$surv)==1) {
	    mat <- cbind(c(fit$n.risk), c(fit$surv), c(fit$std.err*fit$surv))
	    dimnames(mat) <- list(levs, c("n.risk", tname))
	    cat(" Survival at time", fit$time, "\n")
	    surv4.prmatrix(mat)
	    }
	else {
	    for (i in 1:nstrat) {
		cat("       ", levs[i], "\n")
		mat <- cbind(fit$time/scale, fit$n.risk[,i], fit$surv[,i])
		if (!is.null(fit$std.err)) mat<- cbind(mat,
			   fit$std.err[,i] * fit$surv[,i])
		if (!naprint) mat <- mat[!is.na(mat[,3]),,drop=F]
		surv4.prmatrix(mat, rowlab=rep("",nrow(mat)),
				collab=c("Time", "n.risk", tname))
		cat("\n")
		}
	    }
	}
    invisible(fit)
    }
#SCCS 09/27/96 @(#)print.survfit.s	4.8
print.survfit <- function(fit, scale=1, digits = max(options()$digits - 4, 3), ...) {

    if (!is.null(cl<- fit$call)) {
	cat("Call: ")
	dput(cl)
	cat("\n")
	}
    omit <- fit$na.action
    if (length(omit)) cat("  ", naprint(omit), "\n")

    savedig <- options(digits=digits)
    on.exit(options(savedig))
    pfun <- function(stime, surv, n.risk, n.event, lower, upper) {
	#compute the mean, median, se(mean), and ci(median)
	minmin <- function(y, x) {
	     if (any(!is.na(y) & y==.5)) {
	       if (any(!is.na(y) & y <.5))
		 .5*( min(x[!is.na(y) & y==.5]) + min(x[!is.na(y) & y<.5]))
	       else
		 .5*( min(x[!is.na(y) & y==.5]) + max(x[!is.na(y) & y==.5]))
	       }
	     else  min(x[!is.na(y) & y<=.5])
	     }
	n <- length(stime)
	hh <- c(n.event[-n]/(n.risk[-n]*(n.risk[-n]-n.event[-n])), 0)
	nused <- n.risk[1]
	ndead<- sum(n.event)
	dif.time <- c(diff(c(0, stime)), 0)
	if (is.matrix(surv)) {
	    n <- nrow(surv)
	    mean <- dif.time * rbind(1, surv)
	    temp <- (apply(mean[(n+1):2,,drop=F], 2, cumsum))[n:1,,drop=F]
	    varmean <- c(hh %*% temp^2)
	    med <- apply(surv, 2, minmin, stime)
	    if (!is.null(upper)) {
		upper <- apply(upper, 2, minmin, stime)
		lower <- apply(lower, 2, minmin, stime)
		cbind(nused, ndead, apply(mean, 2, sum),
			  sqrt(varmean), med, lower, upper)
		}
	    else {
		cbind(nused, ndead, apply(mean, 2, sum),
			   sqrt(varmean), med)
		}
	    }
	else {
	    mean <- dif.time*c(1, surv)
	    varmean <- sum(rev(cumsum(rev(mean))^2)[-1] * hh)
	    med <- minmin(surv, stime)
	    if (!is.null(upper)) {
		upper <- minmin(upper, stime)
		lower <- minmin(lower, stime)
		c(nused, ndead, sum(mean), sqrt(varmean), med, lower, upper)
		}
	    else {
		c(nused, ndead, sum(mean), sqrt(varmean), med)
		}
	    }
	}

    stime <- fit$time/scale
    surv <- fit$surv
    plab <- c("n", "events", "mean", "se(mean)", "median")
    if (!is.null(fit$conf.int))
	plab2<- paste(fit$conf.int, c("CI", "CI"), sep='')

    #Four cases: strata Y/N  by  ncol(surv)>1 Y/N
    #  Repeat the code, with minor variations, for each one
    if (is.null(fit$strata)) {
	x <- pfun(stime, surv, fit$n.risk, fit$n.event, fit$lower, fit$upper)
	if (is.matrix(x)) {
	    if (is.null(fit$lower)) dimnames(x) <- list(NULL, plab)
	    else                    dimnames(x) <- list(NULL, c(plab, plab2))
	    }
	else {
	    if (is.null(fit$lower)) names(x) <- plab
	    else                    names(x) <- c(plab, plab2)
	    }
	print(x)
	}
    else {   #strata case
	nstrat <- length(fit$strata)
	stemp <- rep(1:nstrat,fit$strata)
	x <- NULL
	for (i in unique(stemp)) {
	    who <- (stemp==i)
	    if (is.matrix(surv)) {
		temp <- pfun(stime[who], surv[who,,drop=F],
			  fit$n.risk[who], fit$n.event[who],
			  fit$lower[who,,drop=F], fit$upper[who,,drop=F])
		x <- rbind(x, temp)
		}
	    else  {
		temp <- pfun(stime[who], surv[who], fit$n.risk[who],
			  fit$n.event[who], fit$lower[who], fit$upper[who])
		x <- rbind(x, temp)
		}
	    }
	temp <- names(fit$strata)
	if (nrow(x) > length(temp)) {
	    nrep <- nrow(x)/length(temp)
	    temp <- rep(temp, rep(nrep, length(temp)))
	    }
	if (is.null(fit$lower)) dimnames(x) <- list(temp, plab)
	else                    dimnames(x) <- list(temp, c(plab, plab2))
	print(x)
	}
    invisible(fit)
    }
#SCCS @(#)print.survreg.s	4.8 11/19/92
print.survreg <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        }
    if (!is.null(x$fail)) {
	cat(" Survreg failed.", x$fail, "\n")
	return(invisible(x))
	}
    coef <- c(x$coef, x$parms[!x$fixed])
    if(any(nas <- is.na(coef))) {
	if(is.null(names(coef))) names(coef) <- paste("b", 1:length(coef), sep = "")
        cat("\nCoefficients: (", sum(nas), 
            " not defined because of singularities)\n", sep = "")
        }
    else cat("\nCoefficients:\n")
    print(coef, ...)
    rank <- x$rank
    if(is.null(rank))
        rank <- sum(!nas)
    nobs <- length(x$residuals)
    rdf <- x$df.resid
    if(is.null(rdf))
        rdf <- nobs - rank
    omit <- x$na.action
    if (length(omit))
	cat("n=", nobs, " (", naprint(omit), ")\n", sep="")
##	sd <- survreg.distributions[[x$family[1]]]
    sd <- survreg.distributions[[x$family[1]$name]]
    cat("\n", sd$print(x$parms, x$fixed), "\n", sep='')
    cat("Degrees of Freedom:", nobs, "Total;", rdf, "Residual\n")
    cat("Residual Deviance:", format(x$deviance), "\n")
    invisible(x)
    }
#SCCS  @(#)pyears.s	4.4 07/29/96
pyears <- function(formula=formula(data), data=sys.frame(sys.parent()),
	weights, subset, na.action,
	ratetable=survexp.us, scale=365.25,  expected=c('event', 'pyears'),
	model=F, x=F, y=F) {

###<TSL>    expect <- match.arg(expect)
    expect <- match.arg(expected)
    call <- match.call()
    m <- match.call(expand=F)
    m$ratetable <- m$model <- m$x <- m$y <- m$scale<- m$expected<-NULL
###<TSL>
###    Terms <- if(missing(data)) terms(formula, 'ratetable')
###	     else              terms(formula, 'ratetable',data=data)
    Terms <-terms(formula, 'ratetable')
### </TSL>
    if (any(attr(Terms, 'order') >1))
	    stop("Pyears cannot have interaction terms")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.frame(sys.parent()))

    Y <- model.extract(m, 'response')
    if (is.null(Y)) stop ("Follow-up time must appear in the formula")
    if (!is.Surv(Y)){
	if (any(Y <0)) stop ("Negative follow up time")
	Y <- as.matrix(Y)
	if (ncol(Y) >2) stop("Y has too many columns")
	if (ncol(Y)==2 && any(Y[,2] <= Y[,1]))
	    stop("stop time must be > start time")
	}
    n <- nrow(Y)

    weights <- model.extract(m, 'weights')

    rate <- attr(Terms, "specials")$ratetable
    if (length(rate) >1 )
	stop ("Can have only 1 ratetable() call in a formula")
    else if (length(rate)==1) {
	ovars <- (dimnames(attr(Terms, 'factors'))[[1]])[-c(1, rate)]
	rtemp <- match.ratetable(m[,rate], ratetable)
	R <- rtemp$R
	if (!is.null(rtemp$call)) {  #need to drop some dimensions from ratetable
	    ratetable <- eval(parse(text=temp$call))
	    }
	}
    else {
	ovars <- (dimnames(attr(Terms, 'factors'))[[1]])[-1]
	}

    # Now process the other (non-ratetable) variables
    if (length(ovars)==0)  {
	# no categories!
	X <- rep(1,n)
	ofac <- odim <- odims <- ocut <- 1
	}
    else {
	odim <- length(ovars)
	ocut <- NULL
	odims <- ofac <- double(odim)
	X <- matrix(0, n, odim)
	outdname <- vector("list", odim)
	for (i in 1:odim) {
	    temp <- m[[ovars[i]]]
	    ctemp <- class(temp)
	    if (!is.null(ctemp) && ctemp=='tcut') {
		X[,i] <- temp
		temp2 <- attr(temp, 'cutpoints')
		odims[i] <- length(temp2) -1
		ocut <- c(ocut, temp2)
		ofac[i] <- 0
		outdname[[i]] <- attr(temp, 'labels')
		}
	    else {
		temp2 <- factor(temp)
		X[,i] <- temp2
		temp3 <- levels(temp2)
		odims[i] <- length(temp3)
		ofac[i] <- 1
		outdname[[i]] <- temp3
		}
	    }
	}

    # Now do the computations
    ocut <-c(ocut,0)   #just in case it were of length 0
    osize <- prod(odims)
    if (length(rate)) {  #include expected
	atts <- attributes(ratetable)
	cuts <- atts$cutpoints
	rfac <- atts$factor
	us.special <- (rfac >1)
	if (any(us.special)) {  #special handling for US pop tables
	    if (sum(us.special) >1)
		stop("Two columns marked for special handling as a US rate table")
	    #slide entry date so that it appears that they were born on Jan 1
	    cols <- match(c("age", "year"), atts$dimid)
	    if (any(is.na(cols))) stop("Ratetable does not have expected shape")
	    temp <- date.mdy(R[,cols[2]]-R[,cols[1]])
	    R[,cols[2]] <- R[,cols[2]] - mdy.date(temp$month, temp$day, 1960)
	    # Doctor up "cutpoints"
	    temp <- (1:length(rfac))[us.special]
	    nyear <- length(cuts[[temp]])
	    nint <- rfac[temp]       #intervals to interpolate over
	    cuts[[temp]] <- round(approx(nint*(1:nyear), cuts[[temp]],
					    nint:(nint*nyear))$y - .0001)
	    }

	temp <- .C("pyears1",
			as.integer(n),
			as.integer(ncol(Y)),
			as.integer(is.Surv(Y)),
			as.double(Y),
			as.integer(length(atts$dim)),
			as.integer(rfac),
			as.integer(atts$dim),
			as.double(unlist(cuts)),
			ratetable,
			as.double(R),
			as.integer(odim),
			as.integer(ofac),
			as.integer(odims),
			as.double(ocut),
			as.integer(expect=='event'),
			X,
			pyears=double(osize),
			pn    =double(osize),
			pcount=double(if(is.Surv(Y)) osize else 1),
			pexpect=double(osize),
			offtable=double(1))[17:21]
	}
    else {
	temp <- .C('pyears2',
			as.integer(n),
			as.integer(ncol(Y)),
			as.integer(is.Surv(Y)),
			as.double(Y),
			as.integer(odim),
			as.integer(ofac),
			as.integer(odims),
			as.double(ocut),
			X,
			pyears=double(osize),
			pn    =double(osize),
			pcount=double(if(is.Surv(Y)) osize else 1),
			offtable=double(1)) [10:13]
	}

    out <- list(call = call,
		pyears= array(temp$pyears/scale, dim=odims, dimnames=outdname),
		n     = array(temp$pn,     dim=odims, dimnames=outdname),
		offtable = temp$offtable/scale)
    if (length(rate)) {
        out$expected <- array(temp$pexpect, dim=odims, dimnames=outdname)
	if (expect=='pyears') out$expected <- out$expected/scale
	if (!is.null(rtemp$summ)) out$summary <- rtemp$summ
	}
    if (is.Surv(Y))
	out$event    <- array(temp$pcount, dim=odims, dimnames=outdname)
    na.action <- attr(m, "na.action")
    if (length(na.action))  out$na.action <- na.action
    if (model) out$model <- m
    else {
	if (x) out$x <- cbind(X, R)
	if (y) out$y <- Y
	}
    class(out) <- 'pyears'
    out
    }

#SCCS @(#)ratetable.s	4.9 01/31/95
#
# This is a 'specials' function for pyears
#   it is a stripped down version of as.matrix(data.frame(...))
# There is no function to create a ratetable.
# This function has a class, only so that data frame subscripting will work
#
ratetable <- function(...) {
    args <- list(...)
    nargs <- length(args)
    ll <- 1:nargs
    for (i in ll) ll[i] <- length(args[[i]])
    n <- max(ll)
    levlist <- vector("list", nargs)
    x <- matrix(0,n,nargs)
    dimnames(x) <- list(1:n, names(args))
    for (i in 1:nargs) {
	if (ll[i] ==n) {
	    if (!is.numeric(args[[i]])) args[[i]] <- factor(args[[i]])
	    if (is.factor(args[[i]])) {
		levlist[[i]] <- levels(args[[i]])
		x[,i] <- c(args[[i]])
		}
	    else x[,i] <- args[[i]]
	    }
	else if (ll[i] ==1) levlist[i] <- args[i]
	else stop("All arguments to ratetable() must be the same length")
	}
    attr(x, "constants") <- (ll==1) & (n>1)
    attr(x, "levlist")   <- levlist
    attr(x, "class")  <- "ratetable2"
    x
    }

# The two functions below should only be called internally, when missing
#   values cause model.frame to drop some rows
is.na.ratetable2 <- function(x) {
    attributes(x) <- list(dim=dim(x))
    as.vector((1 * is.na(x)) %*% rep(1, ncol(x)) >0)
    }
"[.ratetable2" <- function(x, rows, cols, drop=F) {
    if (!missing(cols)) {
       stop("This should never be called!")
       }
    aa <- attributes(x)
    attributes(x) <- aa[c("dim", "dimnames")]
    y <- x[rows,,drop=F]
    attr(y,'constants') <- aa$constants
    attr(y,'levlist')   <- aa$levlist
    class(y) <- aa$class
    y
    }

#
# Functions to manipulate rate tables
#

"[.ratetable" <- function(x, ..., drop=T) {
    aa <- attributes(x)
    attributes(x) <- aa[c("dim", "dimnames")]
    y <- NextMethod("[",x,  drop=F)
    newdim <- attr(y, 'dim')
    if (is.null(newdim)) stop("Invalid subscript")
    dropped <- (newdim==1)
    if (drop)  change <- (newdim!=aa$dim & !dropped)
    else       change <- (newdim!=aa$dim)

    if (any(change)) {  #dims that got smaller, but not dropped
	newcut <- aa$cutpoints
	for (i in (1:length(change))[change])
	    if (!is.null(newcut[[i]])) newcut[[i]] <-
		(newcut[[i]])[match(dimnames(y)[[i]], aa$dimnames[[i]])]
	aa$cutpoints <- newcut
	}
    if (drop && any(dropped)){
	if (all(dropped)) as.numeric(y)   #single element
	else {
	    #Note that we have to drop the summary function
	    attributes(y) <- list( dim = dim(y)[!dropped],
				   dimnames = dimnames(y)[!dropped],
				   dimid = aa$dimid[!dropped],
				   factor = aa$factor[!dropped],
				   cutpoints =aa$cutpoints[!dropped],
				   class = aa$class)
	    y
	    }
	}
    else {
	aa$dim <- aa$dimnames <- NULL
	attributes(y) <- c(attributes(y), aa)
	y
	}
    }

is.na.ratetable  <- function(x)
    structure(is.na(as.vector(x)), dim=dim(x), dimnames=dimnames(x))

Math.ratetable <- function(x, ...) {
    attributes(x) <- attributes(x)[c("dim", "dimnames")]
    NextMethod(.Generic)
    }

Ops.ratetable <- function(e1, e2) {
    #just treat it as an array
    if (nchar(.Method[1])) attributes(e1) <- attributes(e1)[c("dim","dimnames")]
    if (nchar(.Method[2])) attributes(e2) <- attributes(e2)[c("dim","dimnames")]
    NextMethod(.Generic)
    }

as.matrix.ratetable <- function(x) {
    attributes(x) <- attributes(x)[c("dim", "dimnames")]
    x
    }

print.ratetable <- function(x, ...)  {
    cat ("Rate table with dimension(s):", attr(x, 'dimid'), "\n")
    attributes(x) <- attributes(x)[c("dim", "dimnames")]
### <TSL>
###    NextMethod("print")
    print(x)
### </TSL>
    }
#SCCS 04/14/92 @(#)residuals.coxph.null.s	4.2
residuals.coxph.null <-
  function(object, type=c("martingale", "deviance", "score", "schoenfeld"),
	    ...)
    {
    type <- match.arg(type)
    if (type=='martingale' || type=='deviance') NextMethod("residuals",object)
    else stop(paste("\'", type, "\' residuals are not defined for a null model",
			sep=""))
    }
# SCCS @(#)residuals.coxph.s	4.26 11/09/95
residuals.coxph <-
  function(object, type=c("martingale", "deviance", "score", "schoenfeld",
			  "dfbeta", "dfbetas", "scaledsch"),
	    collapse=F, weighted=F)
    {
###    type <- match.arg(type)
    type <- match.arg(type,c("martingale", "deviance", "score", "schoenfeld","dfbeta", "dfbetas", "scaledsch"))
    otype <- type
    if (type=='dfbeta' || type=='dfbetas') type <- 'score'
    if (type=='scaledsch') type<-'schoenfeld'
    n <- length(object$residuals)
    rr <- object$residual
    y <- object$y
    x <- object$x
    vv <- object$naive.var
    if (is.null(vv)) vv <- object$var
    weights <- object$weights
    strat <- object$strata
    method <- object$method
    if (method=='exact' && (type=='score' || type=='schoenfeld'))
	stop(paste(type, 'residuals are not available for the exact method'))

    if (type == 'martingale') rr <- object$residual
    else {
	# I need Y, and perhaps the X matrix (and strata)
	Terms <- object$terms
	if (!inherits(Terms, 'terms'))
		stop("invalid terms component of object")
	strats <- attr(Terms, "specials")$strata
	if (is.null(y)  ||  (is.null(x) && type!= 'deviance')) {
	    temp <- coxph.getdata(object, y=T, x=T, strata=T)
	    y <- temp$y
	    x <- temp$x
	    if (length(strats)!=0) strat <- temp$strata
	    }

	ny <- ncol(y)
	status <- y[,ny,drop=T]

	if (type != 'deviance') {
            ##as.numeric() to codes()
	    nstrat <- if (is.factor(strat)) codes(strat) else as.numeric(strat)
	    nvar <- ncol(x)
	    if (is.null(strat)) {
		ord <- order(y[,ny-1], -status)
		newstrat <- rep(0,n)
		}
	    else {
		ord <- order(nstrat, y[,ny-1], -status)
		newstrat <- c(diff(as.numeric(nstrat[ord]))!=0 ,1)
		}
	    newstrat[n] <- 1

	    # sort the data
	    x <- x[ord,]
	    y <- y[ord,]
	    score <- exp(object$linear.predictor)[ord]
	    if (is.null(weights)) {weights <- rep(1,n); weighted <- F}
	    else                  weights <- weights[ord]
	    }
	}

    #
    # Now I have gotton the data that I need-- do the work
    #
    if (type=='schoenfeld') {
	if (ny==2)  y <- cbind(-1,y)
	temp <- .C("coxscho", n=as.integer(n),
			    as.integer(nvar),
			    as.double(y),
			    resid= x,
			    score * weights,
			    as.integer(newstrat),
			    as.integer(method=='efron'),
			    double(3*nvar))

	deaths <- y[,3]==1

	if (nvar==1) rr <- temp$resid[deaths]
	else rr <- matrix(temp$resid[deaths,], ncol=nvar) #pick rows, and kill attr
	if (length(strats)!=0) attr(rr, "strata")  <- table((strat[ord])[deaths])
	time <- c(y[deaths,2])  # 'c' kills all of the attributes
	if (is.matrix(rr)) dimnames(rr)<- list(time, names(object$coef))
	else               names(rr) <- time

	if (otype=='scaledsch') {
	    ndead <- sum(deaths)
	    coef <- ifelse(is.na(object$coef), 0, object$coef)
	    if (nvar==1) rr <- rr*as.vector(vv) *ndead + coef
	    else         rr <- rr %*%vv * ndead +
						outer(rep(1,nrow(rr)),coef)
	    }
	return(rr)
	}

    if (type=='score') {
	if (ny==2) {
	    resid <- .C("coxscore", as.integer(n),
				as.integer(nvar),
				as.double(y),
				x=x,
				as.integer(newstrat),
				score,
				weights,
				as.integer(method=='efron'),
				resid= double(n*nvar),
				double(2*nvar))$resid
	    }
	else {
	    resid<- .C("agscore",
				as.integer(n),
				as.integer(nvar),
				as.double(y),
				x,
				as.integer(newstrat),
				score,
				weights,
				as.integer(method=='efron'),
				resid=double(n*nvar),
				double(nvar*6))$resid
	    }
	if (nvar >1) {
	    rr <- matrix(0, n, nvar)
	    rr[ord,] <- matrix(resid, ncol=nvar)
	    dimnames(rr) <- list(names(object$resid), names(object$coef))
	    }
	else rr[ord] <- resid

	if      (otype=='dfbeta') {
	    if (is.matrix(rr)) rr <- rr %*% vv
	    else               rr <- rr * vv
	    }
	else if (otype=='dfbetas') {
	    if (is.matrix(rr))  rr <- (rr %*% vv) %*% diag(sqrt(1/diag(vv)))
	    else                rr <- rr * sqrt(vv)
	    }
	}

    #
    # Multiply up by case weights, if requested
    #
    if (!is.null(weights) & weighted) {
	weights[ord] <- weights
	rr <- rr * weights
	}

    #Expand out the missing values in the result
    if (!is.null(object$na.action)) {
	rr <- naresid(object$na.action, rr)
	if (is.matrix(rr)) n <- nrow(rr)
	else               n <- length(rr)
	if (type=='deviance') status <- naresid(object$na.action, status)
	}

    # Collapse if desired
    if (!missing(collapse)) {
	if (length(collapse) !=n) stop("Wrong length for 'collapse'")
	rr <- rowsum(rr, collapse)
	}

    # Deviance residuals are computed after collapsing occurs
    if (type=='deviance')
	sign(rr) *sqrt(-2* (rr+
			      ifelse(status==0, 0, status*log(status-rr))))
    else rr
    }
residuals.survreg <-
function(object, type = c("deviance", "pearson", "working", "matrix"))
{
    type <- match.arg(type)
    rr <- switch(type,
	    working = object$residuals,
	    pearson = sqrt(object$weights) * object$residuals,
	    deviance = object$dresiduals,
	    matrix= {
		eta <- object$linear.predictors
		n   <- length(eta)
		y   <- object$y
		if (is.null(y))
		    stop("Program deficiency-- must have y=T for matrix residuals")
		dist<- match(object$dist, c("extreme", "logistic",
					    "gaussian", "cauchy"))
		temp <-.C("survreg_g", as.integer(n),
				as.double(y),
				as.integer(ncol(y)),
				eta,
				as.double(object$parms),
				deriv=matrix(double(n*6), ncol=6),
				as.integer(6),
				as.integer(dist))$deriv
		dimnames(temp) <- list(names(object$residuals),
				       c("loglik", "eta'", "eta''", "sig'",
					 "sig''", "eta'sig'"))
		temp
		}
	    )

    #Expand out the missing values in the result
    if (!is.null(object$na.action))
	 naresid(object$na.action, rr)
    else rr
    }

# @(#)strata.s	4.9 07/21/93
# Create a strata variable, possibly from many objects
#
strata <- function(..., na.group=F, shortlabel=F) {
    words <- as.character((match.call())[-1])
    if (!missing(na.group)) words <- words[-1]
    allf <- list(...)
    if(length(allf) == 1 && is.list(ttt <- unclass(allf[[1]]))) {
	    allf <- ttt
	    words <- names(ttt)
	    }
    nterms <- length(allf)
    what <- allf[[1]]
### <TSL>
    if(!is.factor(what))
	    what <- factor(what)
###    levs <- unclass(what) - 1
    levs<-codes(what) -1
### </TSL>
    wlab <- levels(what)
    if (na.group && any(is.na(what))){
	levs[is.na(levs)] <- length(wlab)
	wlab <- c(wlab, "NA")
	}
    if (shortlabel) labs <- wlab
    else            labs <- paste(words[1], wlab, sep='=')
    for (i in (1:nterms)[-1]) {
	what <- allf[[i]]
###<TSL>
	if(!is.factor(what))
		what <- factor(what)
	wlab <- levels(what)
### 	wlev <- unclass(what) - 1
	wlev<-codes(what) -1
### </TSL>
	if (na.group && any(is.na(wlev))){
	    wlev[is.na(wlev)] <- length(wlab)
	    wlab <- c(wlab, "NA")
	    }
	if (!shortlabel) wlab <- format(paste(words[i], wlab, sep='='))
	levs <- wlev + levs*(length(wlab))
	labs <- paste(rep(labs, rep(length(wlab), length(labs))),
		      rep(wlab, length(labs)), sep=', ')
	}
    levs <- levs + 1
    ulevs <- sort(unique(levs[!is.na(levs)]))
    levs <- match(levs, ulevs)
    labs <- labs[ulevs]
### <TSL>	class(levs) <- "factor"
	levs<-as.factor(levs)
###</TSL>
    levels(levs) <- labs
    levs
    }
#SCCS 09/27/96 @(#)summary.coxph.s	4.6
summary.coxph <-
 function(cox, table = T, coef = T, conf.int = 0.95, scale = 1,
			digits = max(options()$digits - 4, 3))
    {
    if (!is.null(cl<- cox$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
    if (!is.null(cox$fail)) {
	cat(" Coxreg failed.", cox$fail, "\n")
	return()
	}
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    omit <- cox$na.action
    if (length(omit))
	cat("  n=", cox$n, " (", naprint(omit), ")\n", sep="")
    else cat("  n=", cox$n, "\n")
    if (length(cox$icc))
	cat("  robust variance based on", cox$icc[1],
	    "groups, intra-class correlation =", format(cox$icc[2:3]), "\n")
    if (is.null(cox$coef)) {   # Null model
	cat ("   Null model\n")
	return()
	}

    beta <- cox$coef
    nabeta <- !(is.na(beta))          #non-missing coefs
    beta2 <- beta[nabeta]
    if(is.null(beta) | is.null(cox$var))
        stop("Input is not valid")

    if (is.null(cox$naive.var)) {
	se <- sqrt(diag(cox$var))
	wald.test <-  sum(beta2 * solve(cox$var[nabeta,nabeta], beta2))
	}
    else {
	nse <- sqrt(diag(cox$naive.var))        #naive se
	se <- sqrt(diag(cox$var))
	wald.test <-  sum(beta2 * solve(cox$var[nabeta,nabeta], beta2))
	}
    if(coef) {
	if (is.null(cox$naive.var)) {
	    tmp <- cbind(beta, exp(beta), se, beta/se,
		   signif(1 - pchisq((beta/ se)^2, 1), digits -1))
	    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
		"se(coef)", "z", "p"))
	    }
	else {
	    tmp <- cbind(beta, exp(beta), nse, se, beta/se,
		   signif(1 - pchisq((beta/ se)^2, 1), digits -1))
	    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
		"se(coef)", "robust se", "z", "p"))
	    }
        cat("\n")
        surv4.prmatrix(tmp)
        }
    if(conf.int) {
        z <- qnorm((1 + conf.int)/2, 0, 1)
        beta <- beta * scale
        se <- se * scale
        tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se),
            exp(beta + z * se))
        dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)",
            paste("lower .", round(100 * conf.int, 2), sep = ""),
            paste("upper .", round(100 * conf.int, 2), sep = "")))
        cat("\n")
        surv4.prmatrix(tmp)
        }
    logtest <- -2 * (cox$loglik[1] - cox$loglik[2])
    sctest <- cox$score
    df <- length(beta2)
    cat("\n")
    cat("Rsquare=", format(round(1-exp(-logtest/cox$n),3)),
	"  (max possible=", format(round(1-exp(2*cox$loglik[1]/cox$n),3)),
	")\n" )
    cat("Likelihood ratio test= ", format(round(logtest, 2)), "  on ",
	df, " df,", "   p=", format(1 - pchisq(logtest, df)),
	"\n", sep = "")
    cat("Wald test            = ", format(round(wald.test, 2)), "  on ",
	df, " df,", "   p=", format(1 - pchisq(wald.test, df)),
	"\n", sep = "")
    cat("Efficient score test = ", format(round(sctest, 2)), "  on ", df,
        " df,", "   p=", format(1 - pchisq(sctest, df)), "\n\n", sep = 
        "")
    if (!is.null(cox$naive.var))
      cat("   (Note: the likelihood ratio and efficient score tests",
	  " assume independence of the observations).\n")
    invisible()
    }
#SCCS 02/28/95 @(#)summary.survfit.s	1.7
summary.survfit <- function(fit, times, censored=F, scale=1, ...) {
    if (!inherits(fit, 'survfit'))
	    stop("Invalid data")

    n <- length(fit$time)
    stime <- fit$time/scale
    if (is.null(fit$strata)) {
	stemp <- rep(1,n)
	nstrat <- 1
	}
    else {
	nstrat <- length(fit$strata)
	stemp <- rep(1:nstrat,fit$strata)
	}

    surv <- as.matrix(fit$surv)
    if (is.null(fit$std.err)) std.err <- NULL
    else                      std.err <- fit$std.err * surv

    if (!is.null(fit$lower)) {
	lower <- as.matrix(fit$lower)
	upper <- as.matrix(fit$upper)
	}

    if (missing(times)) {
	if (censored) {
	    times <- stime
	    n.risk<- fit$n.risk
	    n.event <- fit$n.event
	    }
	else {
	    who    <- (fit$n.event >0)
	    times  <-  stime[who]
	    n.risk <-  fit$n.risk[who]
	    n.event <- fit$n.event[who]
	    stemp <- stemp[who]
	    surv <- surv[who,,drop=F]
	    if (!is.null(std.err)) std.err <- std.err[who,,drop=F]
	    if (!is.null(fit$lower)) {
		lower <- lower[who,,drop=F]
		upper <- upper[who,,drop=F]
		}
	    }
	}

    else {  #this case is much harder
	if (any(times<0)) stop("Invalid time point requested")
        if(max(fit$time) < min(times))
            stop("Requested times are all beyond the end of the survival curve")
	if (length(times) >1 )
	    if (any(diff(times)<0)) stop("Times must be in increasing order")

	temp <- .C("survindex2", as.integer(n),
				  as.double(stime),
				  as.integer(stemp),
				  as.integer(length(times)),
				  as.double(times),
				  as.integer(nstrat),
				  indx = integer(nstrat*length(times)),
				  indx2= integer(nstrat*length(times)) )
	keep <- temp$indx >=0
	indx <- temp$indx[keep]
	ones <- (temp$indx2==1)[keep]
	ties <- (temp$indx2==2)[keep]  #data set time === requested time

	times <- rep(times, nstrat)[keep]
	n.risk <- fit$n.risk[indx+1 - (ties+ones)]
	surv   <- surv[indx,,drop=F];   surv[ones,] <- 1
	if (!is.null(std.err)) {
	    std.err<- std.err[indx,,drop=F]
	    std.err[ones,] <-0
	    }
	fit$n.event[stime>max(times)] <- 0
	n.event <- (cumsum(c(0,fit$n.event)))[ifelse(ones, indx, indx+1)]
	n.event<-  diff(c(0, n.event))

	if (!is.null(fit$lower)) {
	    lower <- lower[indx,,drop=F];  lower[ones,] <- 1;
	    upper <- upper[indx,,drop=F];  upper[ones,] <- 1;
	    }

	stemp <- stemp[indx]
	}

    ncurve <- ncol(surv)
    temp <- list(surv=surv, time=times, n.risk=n.risk, n.event=n.event,
			conf.int=fit$conf.int)
    if (ncurve==1) {
	temp$surv <- drop(temp$surv)
	if (!is.null(std.err)) temp$std.err <- drop(std.err)
	if (!is.null(fit$lower)) {
	    temp$lower <- drop(lower)
	    temp$upper <- drop(upper)
	    }
	}
    else {
	if (!is.null(std.err)) temp$std.err <- std.err
	if (!is.null(fit$lower)) {
	    temp$lower <- lower
	    temp$upper <- upper
	    }
	}
    if (!is.null(fit$strata))
	temp$strata <- factor(stemp,
	    labels = names(fit$strata)[sort(unique(stemp))])
    temp$call <- fit$call
    if (!is.null(fit$na.action)) temp$na.action <- fit$na.action
    class(temp) <- 'summary.survfit'
    temp
    }
# SCCS @(#)summary.survreg.s	4.9  12/30/92
summary.survreg<- function(object, correlation = T)
{
    if (!is.null(object$fail)) {
	warning(" Survreg failed.", x$fail, "   No summary provided\n")
	return(invisible(object))
	}
    wt <- object$weights
    fparms <- object$fixed
    coef <- c(object$coef, object$parms[!fparms])
    resid <- object$residuals
    dresid <- object$dresiduals
    n <- length(resid)
    p <- sum(!is.na(coef))
    if(p==0) {
        warning("This model has zero rank --- no summary is provided")
        return(object)
        }
    nsingular <- length(coef) - p
    rdf <- object$df.resid
    if(is.null(rdf))
        rdf <- n - p
# This doesn't seem to do anything
#    R <- object$R   #check for rank deficiencies
#    if(p < max(dim(R)))
#        R <- R[1:p,     #coded by pivoting
#        1:p]
    if(!is.null(wt)) {
        wt <- wt^0.5
        resid <- resid * wt
        excl <- wt == 0
        if(any(excl)) {
            warning(paste(sum(excl), 
                "rows with zero weights not counted"))
            resid <- resid[!excl]
            if(is.null(object$df.residual))
                rdf <- rdf - sum(excl)
            }
        }
    famname <- unlist(object$family["name"])
    if(is.null(famname))
	famname <- "gaussian"
    nas <- is.na(coef)
    cnames <- names(coef[!nas])
    coef <- matrix(rep(coef[!nas], 4), ncol = 4)
    dimnames(coef) <- list(cnames, c("Value", "Std. Error", "z value", "p"))
    stds <- sqrt(diag(object$var[!nas,!nas,drop=F]))
    coef[, 2] <- stds
    coef[, 3] <- coef[, 1]/stds
    coef[, 4] <- 2*pnorm(-abs(coef[,3]))
    if(correlation && sum(!nas)>1 ) {
	correl <- diag(1/stds) %*% object$var[!nas, !nas] %*% diag(1/stds)
        dimnames(correl) <- list(cnames, cnames)
        }
    else correl <- NULL
    ocall <- object$call
    if(!is.null(form <- object$formula)) {
        if(is.null(ocall$formula))
	    ocall <- match.call(get("survreg"), ocall)
        ocall$formula <- form
        }
    sd <- survreg.distributions[[famname]]
    pprint<- paste(sd$name, 'distribution:', sd$print(object$parms, fparms))
    structure(list(call = ocall, terms = object$terms, coefficients = coef,
	scale= NULL, df = c(p, rdf), deviance.resid = dresid,
	var=object$var, correlation = correl, deviance = deviance(object),
	null.deviance = object$null.deviance, iter = object$iter,
	nas = nas, parms=pprint, loglik=object$loglik[2]),
	class = "summary.survreg")
    }
#SCCS @(#)survdiff.fit.s	1.1 01/07/96
survdiff.fit <- function(y, x, strat, rho=0) {
    #
    # This routine is almost always called from survdiff
    #  If called directly, remember that it does no error checking
    #
    n <- length(x)
    if (ncol(y) !=2) stop ("Invalid y matrix")
    if (nrow(y) !=n | length(x) !=n) stop("Data length mismatch")

    ngroup <- length(unique(x))
    if (ngroup <2) stop ("There is only 1 group")
    if (is.category(x)) x <- codes(x)
    else x <- match(x, unique(x))

    if (missing(strat)) strat <- rep(1,n)
    else strat <- as.numeric(as.factor(strat))
    nstrat <- length(unique(strat))
    if (length(strat) !=n) stop("Data length mismatch")

    ord <- order(strat, y[,1], -y[,2])
    strat2 <- c(1*(diff(strat[ord])!=0), 1)

    xx <- .C("survdiff2", as.integer(n),
		   as.integer(ngroup),
		   as.integer(nstrat),
		   as.double(rho),
		   as.double(y[ord,1]),
		   as.integer(y[ord,2]),
		   as.integer(x[ord]),
		   as.integer(strat2),
		   observed = double(ngroup*nstrat),
		   expected = double(ngroup*nstrat),
		   var.e    = double(ngroup * ngroup),
		   double(ngroup), double(n))

    if (nstrat==1)  list(expected = xx$expected,
			 observed = xx$observed,
			 var      = matrix(xx$var, ngroup, ngroup))
    else            list(expected = matrix(xx$expected, ngroup),
			 observed = matrix(xx$observed, ngroup),
			 var      = matrix(xx$var, ngroup, ngroup))
    }
#SCCS 01/07/96 @(#)survdiff.s	4.10
survdiff <- function(formula, data, subset, na.action, rho=0) {
    call <- match.call()
    m <- match.call(expand=F)
    m$rho <- NULL
###
###    Terms <- if(missing(data)) terms(formula, 'strata')
###	     else              terms(formula, 'strata', data=data)
    Terms <-terms(formula, 'strata')
### </TSL>
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.frame(sys.parent()))
    y <- model.extract(m, response)
    if (!inherits(y, "Surv")) stop("Response must be a survival object")
    if (attr(y, 'type') != 'right') stop("Right censored data only")
    ny <- ncol(y)
    n <- nrow(y)

    offset<- attr(Terms, "offset")
    if (!is.null(offset)) {
	#one sample test
	offset <- as.numeric(m[[offset]])
###<TSL>
###	if (length(Terms)>0) stop("Cannot have both an offset and groups")
	if (length(attr(Terms,"term.labels"))>0) stop("Cannot have both an offset and groups")
###</TSL>
	if (any(offset <0 | offset >1))
	    stop("The offset must be a survival probability")
	expected <- sum(1-offset)
	observed <- sum(y[,ny])
	if (rho!=0) {
	    num <- sum(1/rho - ((1/rho + y[,ny])*offset^rho))
	    var <- sum(1- offset^(2*rho))/(2*rho)
	    }
	else {
	    var <-  sum(-log(offset))
	    num <-  var - observed
	    }
	chi <- num*num/var
	rval <-list(n= n, obs = observed, exp=expected, var=var,
			chisq= chi)
	}

    else { #k sample test
	strats <- attr(Terms, "specials")$strata
	if (length(strats)) {
	    temp <- untangle.specials(Terms, 'strata', 1)
	    dropx <- temp$terms
	    if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
	    else strata.keep <- strata(m[,temp$vars], shortlabel=T)
	    }
	else strata.keep <- rep(1,nrow(m))

	#Now create the group variable
	if (length(strats)) ll <- attr(Terms[-dropx], 'term.labels')
	else                ll <- attr(Terms, 'term.labels')
	if (length(ll) == 0) stop("No groups to test")
	else groups <- strata(m[ll])

	fit <- survdiff.fit(y, groups, strata.keep, rho)
	if (is.matrix(fit$observed)){
	    otmp <- apply(fit$observed,1,sum)
	    etmp <- apply(fit$expected,1,sum)
	    }
	else {
	    otmp <- fit$observed
	    etmp <- fit$expected
	    }
	df   <- (etmp >0)                #remove groups with exp=0
	if (sum(df) <2) chi <- 0         # No test, actually
	else {
	    temp2 <- ((otmp - etmp)[df])[-1]
	    vv <- (fit$var[df,df])[-1,-1, drop=F]
	    chi <- sum(solve(vv, temp2) * temp2)
	    }

	rval <-list(n= table(groups), obs = fit$observed,
		    exp = fit$expected, var=fit$var,  chisq=chi)
	if (length(strats)) rval$strata <- table(strata.keep)
	}

    na.action <- attr(m, "na.action")
    if (length(na.action)) rval$na.action <- na.action
    rval$call <- call
    attr(rval, "class") <- 'survdiff'
    rval
    }
# SCCS @(#)survexp.cfit.s	4.5 03/14/95
#
#  Do expected survival based on a Cox model
#   A fair bit of the setup work is identical to survfit.coxph, i.e.,
#     to reconstruct the data frame
#
#  The execution path for individual survival is completely separate, and
#    a whole lot simpler.
#
survexp.cfit <- function(x, y, death, individual, cox, se.fit, method) {
    if (!is.matrix(x)) stop("x must be a matrix")

    #
    # If it is individual survival, things are fairly easy
    #    (the parent routine has guarranteed NO strata in the Cox model
    #
    if (individual) {
	fit <- survfit.coxph(cox, se.fit=F)
	risk <- x[,-1,drop=F] %*% cox$coef  -  sum(cox$coef *cox$means)
	nt <- length(fit$time)
	surv <- approx(-c(0,fit$time), c(1,fit$surv), -y,
				method='constant', rule=2, f=1)$y
	return(list(times=y, surv=c(surv^(exp(risk)))))
	}

    # Otherwise, get on with the real work
    temp <- coxph.getdata(cox, y=T, x=se.fit, strata=F)
    cy <- temp$y
    cx <- temp$x
    cn <- nrow(cy)
    nvar <- length(cox$coef)

    if (ncol(x) != (1+ nvar))
	stop("x matrix does not match the cox fit")

    ngrp <- max(x[,1])
    if (!is.logical(death)) stop("Invalid value for death indicator")

    if (missing(method))
	method <- (1 + 1*(cox$method=='breslow') +2*(cox$method=='efron')
		     + 10*(death))
    else stop("Program unfinished")

    #
    # Data appears ok so proceed
    #  First sort the old data set
    # Also, expand y to (start, stop] form.  This leads to slower processing,
    #  but I only have to program one case instead of 2.
    if (ncol(cy) ==2) cy <- cbind(-1, cy)
    ord <- order(cy[,2], -cy[,3])
    cy  <- cy[ord,]
    score <- exp(cox$linear.predictors[ord])
    if (se.fit) cx <- cx[ord,]
    else  cx <- 0   #dummy, for .C call


    #
    # Process the new data
    #
    if (missing(y) || is.null(y)) y <- rep(max(cy[,2]), nrow(x))
    ord <- order(x[,1])
    x[,1] <- x[,1] - min(x[,1])
    n <- nrow(x)
    ncurve <- length(unique(x[,1]))
    npt <- length(unique(cy[cy[,3]==1,2]))  #unique death times
    xxx  <- .C('agsurv3', as.integer(n),
			  as.integer(nvar),
			  as.integer(ncurve),
			  as.integer(npt),
			  as.integer(se.fit),
			  as.double(score),
			  y = as.double(y[ord]),
			  x[ord,],
			  cox$coef,
			  cox$var,
			  cox$means,
			  as.integer(cn),
			  cy = cy,
			  cx,
			  surv = matrix(0, npt, ncurve),
			  varhaz = matrix(0, npt, ncurve),
			  nrisk  = matrix(0, npt, ncurve),
			  as.integer(method))

    surv <- apply(xxx$surv, 2, cumprod)
    if (se.fit)
	list(surv=surv, n=xxx$nrisk, times=xxx$cy[1:npt],
			se=sqrt(xxx$varhaz)/surv)
    else
	list(surv=surv, n=xxx$nrisk, times=xxx$cy[1:npt,1] )
    }
# SCCS @(#)survexp.fit.s	4.6  05/23/94
#  Actually compute the expected survival for one or more cohorts
#    of subjects.  If each subject is his/her own group, it gives individual
#    survival
survexp.fit <- function(x, y, times, death, ratetable) {
    if (!is.matrix(x)) stop("x must be a matrix")
    if (ncol(x) != (1+length(dim(ratetable))))
	stop("x matrix does not match the rate table")
    atts <- attributes(ratetable)
    rfac <- atts$factor
    if (length(rfac) != ncol(x)-1) stop("Wrong length for rfac")
    ngrp <- max(x[,1])
    times <- sort(unique(times))
    if (any(times <0)) stop("Negative time point requested")
    if (missing(y))  y <- rep(max(times), nrow(x))
    ntime <- length(times)
    if (!is.logical(death)) stop("Invalid value for death indicator")

    cuts <- atts$cutpoints
    us.special <- (rfac >1)
    if (any(us.special)) {  #special handling for US pop tables
	if (sum(us.special) >1)
	    stop("Two columns marked for special handling as a US rate table")
	#slide entry date so that it appears that they were born on Jan 1
	cols <- 1+match(c("age", "year"), attr(ratetable, "dimid"))
	if (any(is.na(cols))) stop("Ratetable does not have expected shape")
	temp <- date.mdy(x[,cols[2]]-x[,cols[1]])
	x[,cols[2]] <- x[,cols[2]] - mdy.date(temp$month, temp$day, 1960)
	# Doctor up "cutpoints"
	temp <- (1:length(rfac))[us.special]
	nyear <- length(cuts[[temp]])
	nint <- rfac[temp]       #intervals to interpolate over
	cuts[[temp]] <- round(approx(nint*(1:nyear), cuts[[temp]],
					nint:(nint*nyear))$y - .0001)
	}
    temp <- .C('pyears3',
		    as.integer(death),
		    as.integer(nrow(x)),
		    as.integer(length(atts$dim)),
		    as.integer(rfac),
		    as.integer(atts$dim),
		    as.double(unlist(cuts)),
		    ratetable,
		    as.double(x),
		    as.double(y),
		    as.integer(ntime),
		    as.integer(ngrp),
		    as.double(times),
		    surv = double(ntime * ngrp),
		    n   = integer(ntime *ngrp))
    if (ntime==1) list(surv=temp$surv, n=temp$n)
    else if (ngrp >1)
	 list(surv=apply(matrix(temp$surv, ntime, ngrp),2,cumprod),
		 n=   matrix(temp$n, ntime, ngrp))
    else list(surv=cumprod(temp$surv), n=temp$n)
    }
#SCCS  @(#)survexp.s	4.21 03/14/95
survexp <- function(formula=formula(data), data=sys.frame(sys.parent()),
	weights, subset, na.action,
	times,  cohort=T,  conditional=F,
	ratetable=survexp.us, scale=1, npoints, se.fit,
	model=F, x=F, y=F) {

    call <- match.call()
    m <- match.call(expand=F)
    m$ratetable <- m$model <- m$x <- m$y <- m$scale<- m$cohort <- NULL
    m$times <- m$conditional <- m$npoints <- m$se.fit <- NULL
###<TSL>
###    Terms <- if(missing(data)) terms(formula, 'ratetable')
###	     else              terms(formula, 'ratetable',data=data)
###    old.as.numeric<-function(x) if (is.factor(x)) codes(x) else as.numeric(x)
 	Terms <- terms(formula, 'ratetable')
###</TSL>
    if (any(attr(Terms, 'order') >1))
	    stop("Survexp cannot have interaction terms")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.frame(sys.parent()))
    n <- nrow(m)

    if (!missing(times)) {
	if (any(times<0)) stop("Invalid time point requested")
	if (length(times) >1 )
	    if (any(diff(times)<0)) stop("Times must be in increasing order")
	}

    Y <- model.extract(m, 'response')
    no.Y <- is.null(Y)
    if (!no.Y) {
	if (is.matrix(Y)) {
	    if (is.Surv(Y) && attr(Y, 'type')=='right') Y <- Y[,1]
	    else stop("Illegal response value")
	    }
	if (any(Y<0)) stop ("Negative follow up time")
	if (missing(npoints)) temp <- unique(Y)
	else                  temp <- seq(min(Y), max(Y), length=npoints)
	if (missing(times)) newtime <- sort(temp)
	else  newtime <- sort(unique(c(times, temp[temp<max(times)])))
	}
    else conditional <- F
    weights <- model.extract(m, 'weights')
    if (!is.null(weights)) warning("Weights ignored")
     rate <- attr(Terms, "specials")$ratetable
    if (length(rate)==0)
	stop("Must have a ratetable() call in the formula")
    if (length(rate) >1 )
	stop ("Can have only 1 ratetable() call in a formula")

    if (no.Y) ovars <- attr(Terms, 'term.labels')[-rate]
    else      ovars <- attr(Terms, 'term.labels')[-(rate-1)]

    if (is.ratetable(ratetable)) {
	israte <- T
	if (no.Y) {
	    if (missing(times))
	       stop("There is no times argument, and no follow-up times are given in the formula")
	    else newtime <- sort(unique(times))
	    Y <- rep(max(times), n)
	    }
	se.fit <- F
	rtemp <- match.ratetable(m[,rate], ratetable)
	R <- rtemp$R
	if (!is.null(rtemp$call)) {  #need to dop some dimensions from ratetable
	    ratetable <- eval(parse(text=rtemp$call))
	    }
       }
    else if (inherits(ratetable, 'coxph')) {
	israte <- F
	Terms <- ratetable$terms
	if (!inherits(Terms, 'terms'))
		stop("invalid terms component of fit")
	if (!is.null(attr(Terms, 'offset')))
	    stop("Cannot deal with models that contain an offset")
	m2 <- data.frame(unclass(m[,rate]))
	strats <- attr(Terms, "specials")$strata
	if (length(strats))
	    stop("survexp cannot handle stratified Cox models")
	R <- model.matrix(delete.response(Terms), m2)[,-1,drop=F]
	if (any(dimnames(R)[[2]] != names(ratetable$coef)))
	    stop("Unable to match new data to old formula")
	if (no.Y) {
	    if (missing(se.fit)) se.fit <- T
	    }
	else se.fit <- F
	}
    else stop("Invalid ratetable argument")

    if (cohort) {
	# Now process the other (non-ratetable) variables
	if (length(ovars)==0)  X <- rep(1,n)  #no categories
	else {
	    odim <- length(ovars)
	    for (i in 1:odim) {
		temp <- m[[ovars[i]]]
		ctemp <- class(temp)
		if (!is.null(ctemp) && ctemp=='tcut')
		    stop("Can't use tcut variables in expected survival")
		}
	    X <- strata(m[ovars])
	    }

	#do the work
	if (israte)
	    temp <- survexp.fit(cbind(as.numeric(X),R), Y, newtime,
			       conditional, ratetable)
	else {
	    temp <- survexp.cfit(cbind(as.numeric(X),R), Y, conditional, F,
			       ratetable, se.fit=se.fit)
	    newtime <- temp$times
	    }
	#package the results
	if (missing(times)) {
	    n.risk <- temp$n
	    surv <- temp$surv
	    if (se.fit) err <- temp$se
	    }
	else {
	    if (israte) keep <- match(times, newtime)
	    else {
		# taken straight out of summary.survfit....
		n <- length(temp$times)
		temp2 <- .C("survindex2", as.integer(n),
					  as.double(temp$times),
					  as.integer(rep(1,n)),
					  as.integer(length(times)),
					  as.double(times),
					  as.integer(1),
					  indx = integer(length(times)),
					  indx2= integer(length(times)) )
		keep <- temp2$indx[temp2$indx>0]
		}

	    if (is.matrix(temp$surv)) {
		surv <- temp$surv[keep,,drop=F]
		n.risk <- temp$n[keep,,drop=F]
		if (se.fit) err <- temp$se[keep,,drop=F]
		}
	    else {
		surv <- temp$surv[keep]
		n.risk <- temp$n[keep]
		if (se.fit) err <- temp$se[keep]
		}
	    newtime <- times
	    }
	newtime <- newtime/scale
	if (length(ovars)) {    #matrix output
	    if (no.Y && israte){ # n's are all the same, so just send a vector
		dimnames(surv) <- list(NULL, levels(X))
		out <- list(call=call, surv=surv, n.risk=c(n.risk[,1]),
			    time=newtime)
		}
	    else {
		#Need a matrix of n's, and a strata component
		out <- list(call=call, surv=surv, n.risk=n.risk,
				time = newtime)
		tstrat <- rep(nrow(surv), ncol(surv))
		names(tstrat) <- levels(X)
		out$strata <- tstrat
		}
	    if (se.fit) out$std.err <- err
	    }
	else {
	     out <- list(call=call, surv=c(surv), n.risk=c(n.risk),
			   time=newtime)
	     if (se.fit) out$std.err <- c(err)
	     }

	na.action <- attr(m, "na.action")
	if (length(na.action))  out$na.action <- na.action
	if (model) out$model <- m
	else {
	    if (x) out$x <- structure(cbind(X, R),
		dimnames=list(row.names(m), c("group", dimid)))
	    if (y) out$y <- Y
	    }
	if (israte && !is.null(rtemp$summ)) out$summ <- rtemp$summ
	if (no.Y) out$method <- 'exact'
	else if (conditional) out$method <- 'conditional'
	else                  out$method <- 'cohort'
	class(out) <- c("survexp", "survfit")
	out
	}

    else { #individual survival
	if (no.Y) stop("For non-cohort, an observation time must be given")
	if (israte)
	    temp <- survexp.fit (cbind(1:n,R), Y, max(Y), T, ratetable)
	else temp<- survexp.cfit(cbind(1:n,R), Y, F, T, ratetable, F)
	xx <- temp$surv
	names(xx) <- row.names(m)
	na.action <- attr(m, "na.action")
	if (length(na.action)) naresid(na.action, xx)
	else xx
	}
    }
#SCCS  @(#)survfit.coxph.null.s	4.10 08/25/94
survfit.coxph.null <-
  function(object, newdata, se.fit=T, conf.int=.95, individual=F,
	    type=c('tsiatis', 'kaplan-meier'),
	    conf.type=c('log', 'log-log', 'plain', 'none'), ...) {
    # May have strata and/or offset terms, linear predictor = offset
    #  newdata doesn't make any sense
    #  This is survfit.coxph with lots of lines removed

    call <- match.call()
    Terms <- terms(object)
    strat <- attr(Terms, "specials")$strata
    n <- object$n
    score <- exp(object$linear.predictor)
    method <- match.arg(type)
    if (!se.fit) conf.type <- 'none'
    else conf.type <- match.arg(conf.type)

    y <- object$y
    stratx <- object$strata
    if (is.null(y) || (length(strat) && is.null(stratx))) {
	# I need the model frame
	m <- model.frame(object)
	if (is.null(stratx)) {
	    temp <- untangle.specials(Terms, 'strata', 1)
	    stratx <- strata(m[temp$vars])
	    }
	if (is.null(y)) y <- model.extract(m, 'response')
	}
    if (is.null(stratx)) stratx <- rep(1,n)
    ny <- ncol(y)
    if (nrow(y) != n) stop ("Mismatched lengths: logic error")

    type <- attr(y, 'type')
    if (type=='counting') {
	ord <- order(stratx, y[,2], -y[,3])
	if (method=='kaplan-meier')
	      stop ("KM method not valid for counting type data")
	}
    else if (type=='right') {
	ord <- order(stratx, y[,1], -y[,2])
	y <- cbind(-1, y)
	}
    else stop("Cannot handle \"", type, "\" type survival data")

    if (length(strat)) {
        ## as.numeric() to codes()
	newstrat <- (codes(stratx))[ord]
	newstrat <- as.integer(c(1*(diff(newstrat)!=0), 1))
	}
    else newstrat <- as.integer(c(rep(0,n-1),1))

    if ( !missing(newdata))
	stop("A newdata argument does not make sense for a null model")

    dimnames(y) <- NULL   #I only use part of Y, so names become invalid
    surv <- .C('agsurv2', as.integer(n),
			  as.integer(0),
			  y = y[ord,],
			  score[ord],
			  strata = newstrat,
			  surv = double(n),
			  varhaz = double(n),
			  double(1),
			  double(0),
			  nsurv = as.integer(method=='kaplan-meier'),
			  double(2),
			  as.integer(1),
			  double(1),
			  newrisk= as.double(1))
    nsurv <- surv$nsurv
    ntime <- 1:nsurv
    tsurv <- surv$surv[ntime]
    tvar  <- surv$varhaz[ntime]
    if (surv$strata[1] <=1)
	temp _ list(time=surv$y[ntime,1],
		 n.risk=surv$y[ntime,2],
		 n.event=surv$y[ntime,3],
		 surv=tsurv )
    else {
	temp <- surv$strata[1:(1+surv$strata[1])]
	tstrat <- diff(c(0, temp[-1])) #n in each strata
	names(tstrat) <- levels(stratx)
	temp _ list(time=surv$y[ntime,1],
		 n.risk=surv$y[ntime,2],
		 n.event=surv$y[ntime,3],
		 surv=tsurv,
		 strata= tstrat)
	}
    if (se.fit) temp$std.err <- sqrt(tvar)

    zval _ qnorm(1- (1-conf.int)/2, 0,1)
    if (conf.type=='plain') {
	temp1 <- temp$surv + zval* temp$std * temp$surv
	temp2 <- temp$surv - zval* temp$std * temp$surv
	temp <- c(temp, list(upper=pmin(temp1,1), lower=pmax(temp2,0),
			conf.type='plain', conf.int=conf.int))
	}
    if (conf.type=='log') {
	xx <- ifelse(temp$surv==0,1,temp$surv)  #avoid some "log(0)" messages
	temp1 <- ifelse(temp$surv==0, 0*temp$std, exp(log(xx) + zval* temp$std))
	temp2 <- ifelse(temp$surv==0, 0*temp$std, exp(log(xx) - zval* temp$std))
	temp <- c(temp, list(upper=pmin(temp1,1), lower=temp2,
			conf.type='log', conf.int=conf.int))
	}
    if (conf.type=='log-log') {
	who <- (temp$surv==0 | temp$surv==1) #special cases
	xx <- ifelse(who, .1,temp$surv)  #avoid some "log(0)" messages
	temp1 <- exp(-exp(log(-log(xx)) + zval*temp$std/log(xx)))
	temp1 <- ifelse(who, temp$surv + 0*temp$std, temp1)
	temp2 <- exp(-exp(log(-log(xx)) - zval*temp$std/log(xx)))
	temp2 <- ifelse(who, temp$surv + 0*temp$std, temp2)
	temp <- c(temp, list(upper=temp1, lower=temp2,
			conf.type='log-log', conf.int=conf.int))
	}

    temp$call <- call
    attr(temp, 'class') <- c("survfit.cox", "survfit")
    temp
    }
#SCCS @(#)survfit.coxph.s	4.15 06/17/93
survfit.coxph <-
  function(object, newdata, se.fit=T, conf.int=.95, individual=F,
	    type=c('tsiatis', 'kaplan-meier'),
	    conf.type=c('log', 'log-log', 'plain', 'none'))
		 {

    if(!is.null((object$call)$weights))
	stop("Survfit cannot (yet) compute the result for a weighted model")
    call <- match.call()
    Terms <- terms(object)
    strat <- attr(Terms, "specials")$strata
### <TSL> Not used
###	    resp <-  attr(Terms, "variables")[attr(Terms, "response")]
    n <- object$n
    nvar <- length(object$coef)
    score <- exp(object$linear.predictor)
    method <- match.arg(type)
    coxmethod <- object$method
    if (!se.fit) conf.type <- 'none'
    else conf.type <- match.arg(conf.type)
    if (length(strat)) temp <- untangle.specials(Terms, 'strata', 1)

    x <- object$x
    y <- object$y
    if (is.null(x) && (length(strat) || se.fit)) {  # I need both X and Y
	stratum <- object$strata
	m <- model.frame(object)
	if (is.null(x)) {   #Both strata and X will be null, or neither
	    if (length(strat)) {
		x <- model.matrix("[.terms"(Terms,-temp$terms), m)[,-1,drop=F]
		stratum <- strata(m[temp$vars])
		}
	    else {
		x <- model.matrix(Terms, m)[,-1,drop=F]
		stratum <- rep(1,n)
		}
	    }
	if (is.null(y)) y <- model.extract(m, 'response')
	}
    else {
	y <- object$y
	if (is.null(y)) {
	    m <- model.frame(object)
	    y <- model.extract(m, 'response')
	    }
	if (length(strat)) stratum <- object$strata
	else               stratum <- factor(rep(1,n))
	}
    if (is.null(x)) x <- matrix(0,n,nvar)

    ny <- ncol(y)
    if (nrow(y) != n) stop ("Mismatched lengths: logic error")
    type <- attr(y, 'type')
    if (type=='counting') {
	ord <- order(stratum, y[,2], -y[,3])
	if (method=='kaplan-meier')
	      stop ("KM method not valid for counting type data")
	}
    else if (type=='right') {
	ord <- order(stratum, y[,1], -y[,2])
	y <- cbind(-1, y)
	}
    else stop("Cannot handle \"", type, "\" type survival data")

    if (length(strat)) {
      ##as.numeric() to codes()
	newstrat <- (codes(stratum))[ord]
	newstrat <- as.integer(c(1*(diff(newstrat)!=0), 1))
	}
    else newstrat <- as.integer(c(rep(0,n-1),1))

    if (individual && !missing(newdata)) stype <- 1
    else {
	stype <- 2
	if (length(strat)) Terms <- "[.terms"(Terms,-temp$terms)  #don't need it
	}
    offset2 <- mean(object$linear.predictors)
    if (!missing(newdata)) {
	m2 <- model.newframe(Terms, newdata, response=(stype==1))
	if (!inherits(m2, 'data.frame'))  {
	    x2 <- as.matrix(m2)
	    if (ncol(x2) != nvar) stop ("Wrong # of variables in new data")
	    n2 <- nrow(x2)
	    if (stype==1) stop("Program error #3")
	    }

	else  {
###	    x2 <- model.matrix(Terms, m2)[,-1,drop=F]
	    if (attr(attr(m2,"terms"),"response"))
	      x2<-model.matrix(Terms,m2)[,-1,drop=F]
	    else
	      x2 <- model.matrix(delete.response(Terms), m2)[,-1,drop=F]
###
	    n2 <- nrow(x2)
	    offset2 <- model.extract(m2, 'offset')
	    if (is.null(offset2)) offset2 <- 0
	    if (stype==1) {
		#
		# The case of an agreg, with a multiple line newdata
		#
		if (length(strat)) {
		    strata2 <- factor(x2[,strat], levels=levels(stratum))
		    x2 <- x2[, -strat, drop=F]
		    }
		else strata2 <- rep(1, nrow(x2))
		y2 <- model.extract(m2, 'response')
		if (attr(y2,'type') != type)
		    stop("Survival type of newdata does not match the fitted model")
		if (nrow(y2) != n2) stop("Wrong # of rows for Y")
		}
	    }
	}
    else x2 <- matrix(object$means, nrow=1)
    n2 <- nrow(x2)
    coef <- ifelse(is.na(object$coef), 0, object$coef)
    newrisk <- exp(c(x2 %*% coef) + offset2 - sum(coef*object$means))

    dimnames(y) <- NULL   #I only use part of Y, so names become invalid
    if (stype==1) {
	surv <- .C("agsurv1", as.integer(n),
			     as.integer(nvar),
			     y[ord,],
			     score,
			     strata=newstrat,
			     surv=double(n*n2),
			     varh=double(n*n2),
			     nsurv=as.integer(2+ 1*(coxmethod=='efron')),
			     x[ord,],
			     double(3*nvar),
			     object$var,
			     y = double(3*n*n2),
			     as.integer(n2),
			     y2,
			     x2,
			     newrisk,
			     as.integer(strata2) )
	ntime <- 1:surv$nsurv
	temp <- (matrix(surv$y, ncol=3))[ntime,]
	temp <- list(time = temp[,1],
		     n.risk= temp[,2],
		     n.event=temp[,3],
		     surv = surv$surv[ntime])
	if (se.fit) temp$std.err <- sqrt(surv$varh[ntime])
	}
    else {
	temp <- ifelse(method=='kaplan-meier', 1,
					2+as.integer(coxmethod=='efron'))
	surv <- .C('agsurv2', as.integer(n),
			      as.integer(nvar* se.fit),
			      y = y[ord,],
			      score[ord],
			      strata = newstrat,
			      surv = double(n*n2),
			      varhaz = double(n*n2),
			      x[ord,],
			      object$var,
			      nsurv = as.integer(temp),
			      double(3*nvar),
			      as.integer(n2),
			      x2,
			      newrisk)
	nsurv <- surv$nsurv
	ntime <- 1:nsurv
	if (n2>1) {
	    tsurv <- matrix(surv$surv[1:(nsurv*n2)], ncol=n2)
	    tvar  <- matrix(surv$varhaz[1:(nsurv*n2)], ncol=n2)
	    dimnames(tsurv) <- list(NULL, dimnames(x2)[[1]])
	    }
	else {
	    tsurv <- surv$surv[ntime]
	    tvar  <- surv$varhaz[ntime]
	    }
	if (surv$strata[1] <=1)
	    temp _ list(time=surv$y[ntime,1],
		     n.risk=surv$y[ntime,2],
		     n.event=surv$y[ntime,3],
		     surv=tsurv )
	else {
	    temp <- surv$strata[1:(1+surv$strata[1])]
	    tstrat <- diff(c(0, temp[-1])) #n in each strata
	    names(tstrat) <- levels(stratum)
	    temp _ list(time=surv$y[ntime,1],
		     n.risk=surv$y[ntime,2],
		     n.event=surv$y[ntime,3],
		     surv=tsurv,
		     strata= tstrat)
	    }
	if (se.fit) temp$std.err <- sqrt(tvar)
	}

    zval _ qnorm(1- (1-conf.int)/2, 0,1)
    if (conf.type=='plain') {
	temp1 <- temp$surv + zval* temp$std * temp$surv
	temp2 <- temp$surv - zval* temp$std * temp$surv
	temp <- c(temp, list(upper=pmin(temp1,1), lower=pmax(temp2,0),
			conf.type='plain', conf.int=conf.int))
	}
    if (conf.type=='log') {
	xx <- ifelse(temp$surv==0,1,temp$surv)  #avoid some "log(0)" messages
	temp1 <- ifelse(temp$surv==0, 0*temp$std, exp(log(xx) + zval* temp$std))
	temp2 <- ifelse(temp$surv==0, 0*temp$std, exp(log(xx) - zval* temp$std))
	temp <- c(temp, list(upper=pmin(temp1,1), lower=temp2,
			conf.type='log', conf.int=conf.int))
	}
    if (conf.type=='log-log') {
	who <- (temp$surv==0 | temp$surv==1) #special cases
	xx <- ifelse(who, .1,temp$surv)  #avoid some "log(0)" messages
	temp1 <- exp(-exp(log(-log(xx)) + zval*temp$std/log(xx)))
	temp1 <- ifelse(who, temp$surv + 0*temp$std, temp1)
	temp2 <- exp(-exp(log(-log(xx)) - zval*temp$std/log(xx)))
	temp2 <- ifelse(who, temp$surv + 0*temp$std, temp2)
	temp <- c(temp, list(upper=temp1, lower=temp2,
			conf.type='log-log', conf.int=conf.int))
	}

    temp$call <- call
    attr(temp, 'class') <- c("survfit.cox", "survfit")
    temp
    }
#SCCS @(#)survfit.km.s	4.12 09/17/96
survfit.km <- function(x, y, casewt=rep(1,n),
	    type=c('kaplan-meier', 'fleming-harrington', 'fh2'),
	    error=c('greenwood', "tsiatis"), se.fit=T,
	    conf.int= .95,
	    conf.type=c('log',  'log-log',  'plain', 'none'),
	    conf.lower=c('usual', 'peto', 'modified'))
    {
    type <- match.arg(type)
    method <- match(type, c("kaplan-meier", "fleming-harrington", "fh2"))

    error <- match.arg(error)
    error.int <- match(error, c("greenwood", "tsiatis"))
    conf.type <- match.arg(conf.type)
    conf.lower<- match.arg(conf.lower)

    ny <- ncol(y)
    n <- nrow(y)

    if (!is.Surv(y)) stop("y must be a Surv object")
    if (!is.factor(x)) stop("x must be a factor")
    if (attr(y, 'type') != 'right') stop("Can only handle right censored data")

    sorted <- (1:n)[order(x, y[,ny-1])]
    y <- y[sorted,]
    ##as.numeric to codes()
    newstrat <- codes(x[sorted])
    newstrat <- as.integer(c(1*(diff(newstrat)!=0), 1))
    if (sum(newstrat) > n/2)
	stop("Number of strata > number of observations/2")
    if (method==3 && any(floor(casewt) != casewt))
	stop("The fh2 method is not valid for fractional case weights")

    storage.mode(y) <- "double"
    dimnames(y) <- NULL
    surv <- .C("survfit2", as.integer(n),
			  y = y,
			  as.integer(ny),
			  as.double(casewt[sorted]),
			  strata= as.integer(newstrat),
			  nstrat= as.integer(method),
			  as.integer(error.int),
			  mark=double(n),
			  surv=double(n),
			  varhaz=double(n),
			  risksum=double(n),
			  ntime = integer(1))
    ntime <- surv$ntime
    if (error.int==1) surv$varhaz[surv$surv==0] <- NA
    ntime <- 1:ntime
    if (surv$nstrat ==1)
	temp _ list(time=surv$y[ntime,1],
		 n.risk=surv$risksum[ntime],
		 n.event=surv$mark[ntime],
		 surv=surv$surv[ntime])
    else {
	temp <- surv$strata[1:surv$nstrat]
	tstrat <- diff(c(0, temp)) #n in each strata
	names(tstrat) <- levels(x)
	temp _ list(time=surv$y[ntime,1],
		 n.risk=surv$risksum[ntime],
		 n.event=surv$mark[ntime],
		 surv=surv$surv[ntime],
		 strata= tstrat)
	}

    if (se.fit) {
	std.err <- sqrt(surv$varhaz[ntime])
	temp$std.err <- std.err
	events <- temp$n.event >0
	n.lag <- rep(c(temp$n.risk[1], temp$n.risk[events]),
	              diff(c(ntime[1], ntime[events], 1+max(ntime))))
	std.low <- switch(conf.lower,
			'usual'   = std.err,
			'peto'    = sqrt((1-temp$surv)/ temp$n.risk),
			'modified'= std.err * sqrt(n.lag/temp$n.risk)
			)
	zval _ qnorm(1- (1-conf.int)/2, 0,1)
	if (conf.type=='plain') {
	    temp1 <- temp$surv + zval* std.err * temp$surv
	    temp2 <- temp$surv - zval* std.low * temp$surv
	    temp <- c(temp, list(upper=pmin(temp1,1), lower=pmax(temp2,0),
			    conf.type='plain', conf.int=conf.int))
	    }
	if (conf.type=='log') {
	    xx <- ifelse(temp$surv==0,1,temp$surv)  #avoid some "log(0)" messages
	    temp1 <- ifelse(temp$surv==0, NA, exp(log(xx) + zval* std.err))
	    temp2 <- ifelse(temp$surv==0, NA, exp(log(xx) - zval* std.low))
	    temp <- c(temp, list(upper=pmin(temp1,1), lower=temp2,
			    conf.type='log', conf.int=conf.int))
	    }
	if (conf.type=='log-log') {
	    who <- (temp$surv==0 | temp$surv==1) #special cases
	    temp3 <- ifelse(temp$surv==0, NA, 1)
	    xx <- ifelse(who, .1,temp$surv)  #avoid some "log(0)" messages
	    temp1 <- exp(-exp(log(-log(xx)) + zval*std.err/log(xx)))
	    temp1 <- ifelse(who, temp3, temp1)
	    temp2 <- exp(-exp(log(-log(xx)) - zval*std.low/log(xx)))
	    temp2 <- ifelse(who, temp3, temp2)
	    temp <- c(temp, list(upper=temp1, lower=temp2,
			    conf.type='log-log', conf.int=conf.int))
	    }
	}
    temp
    }
#SCCS 05/23/95 @(#)survfit.s	4.7
survfit <- function (formular, data, weights, subset, na.action, ...) {
    call <- match.call()
    # Real tricky -- find out if the first arg is "Surv(...)" without
    #  evaluating it.  If this is so, or it is a survival object, turn it
    #  into a formula
###<TSL> even trickier in R
###    if ((mode(call[[2]]) == 'call' &&  call[[2]][[1]] == as.name('Surv')) || inherits(formula, 'Surv'))  
###
	if (((mode(call[[2]])=="call") &&  as.character(call[[2]][[1]])=="Surv") || inherits(formular,'Surv'))
{
	# The dummy function stops an annoying warning message "Looking for
	#  'formula' of mode function, ignored one of mode ..."
	###<TSL> except that it doesn't in R, so I changed the argument name
	xx <- function(x) formula(x)
	formular <- xx(paste(deparse(call[[2]]), 1, sep="~"))
	}

    if (!inherits(formular, 'formula')) {
      ### Rchange
      temp <- UseMethod("survfit")
    }
    else {
	m <- match.call(expand=F)
	m$... <- NULL

	Terms <- terms(formular, 'strata')
	ord <- attr(Terms, 'order')
	if (length(ord) & any(ord !=1))
	    stop("Interaction terms are not valid for this function")
	m$formula <- Terms
        m$formular<-NULL
	m[[1]] <- as.name("model.frame")
	m <- eval(m, sys.frame(sys.parent()))

	n <- nrow(m)
	Y <- model.extract(m, response)
	casewt <- model.extract(m, "weights")
	# The second line below works around a bug in Splus 3.0.1, which later
	#    went away, i.e., casewt is returned as an unevaluated arg.
	if (is.null(casewt)) casewt <- rep(1,n)
	else if (mode(casewt)=='argument') casewt <- eval(casewt[[1]])

	if (!is.null(attr(Terms, 'offset'))) warning("Offset term ignored")

	ll <- attr(Terms, 'term.labels')
	if (length(ll) == 0) X <- factor(rep(1,n))
	else {
	    temp <-  rep(1, length(ll))
	    strat <- untangle.specials(Terms, 'strata',1)$terms
	    if (length(strat)) temp[strat] <- 0
	    lname <- ifelse(temp==1, paste(ll,'=', sep=''), "")
	    temp <- paste("'",lname, "', ","as.character(m[['", ll, "']])", sep='')
	    temp <- paste(temp, collapse=", ', ',")
	    temp <- paste("paste(", temp, ",sep='')")
## Rchange added [[1]]
	    X <- factor(eval(parse(text=temp)[[1]]))
	    }

	temp <- survfit.km(X, Y, casewt, ...)
	attr(temp, "class") <- "survfit"
	if (!is.null(attr(m, 'na.action'))) temp$na.action <- attr(m, 'na.action')
	}
    temp$call <- call
    temp
    }

# The subscript function is bundled in here, although used most
#  often in plotting

"[.survfit" <- function(fit,i,j, drop=F) {
    if (is.null(fit$strata)) {
	if (is.matrix(fit$surv)) {
	    fit$surv <- fit$surv[,i,drop=drop]
	    if (!is.null(fit$std.err)) fit$std.err <- fit$std.err[,i,drop=drop]
	    if (!is.null(fit$upper)) fit$upper <- fit$upper[,i,drop=drop]
	    if (!is.null(fit$lower)) fit$lower <- fit$lower[,i,drop=drop]
	    }
	else warning("Survfit object has only a single survival curve")
	}
    else {
	if (missing(i)) keep <- seq(along=fit$time)
	else {
	    if (is.character(i)) strat <- rep(names(fit$strata), fit$strata)
	    else                 strat <- rep(1:length(fit$strata), fit$strata)
	    keep <- seq(along=strat)[match(strat, i, nomatch=0)>0]
	    fit$strata  <- fit$strata[i]
	    fit$time    <- fit$time[keep]
	    fit$n.risk  <- fit$n.risk[keep]
	    fit$n.event <- fit$n.event[keep]
	    }
	if (is.matrix(fit$surv)) {
	    if (missing(j)) {
		fit$surv <- fit$surv[keep,drop=drop]
		if (!is.null(fit$std.err)) fit$std.err <- fit$std.err[keep,drop=drop]
		if (!is.null(fit$upper)) fit$upper <- fit$upper[keep,drop=drop]
		if (!is.null(fit$lower)) fit$lower <- fit$lower[keep,drop=drop]
		}
	    else {
		fit$surv <- fit$surv[keep,j]
		if (!is.null(fit$std.err)) fit$std.err <- fit$std.err[keep,j]
		if (!is.null(fit$upper)) fit$upper <- fit$upper[keep,j]
		if (!is.null(fit$lower)) fit$lower <- fit$lower[keep,j]
		}
	    }
	else {
	    fit$surv <- fit$surv[keep]
	    if (!is.null(fit$std.err)) fit$std.err <- fit$std.err[keep]
	    if (!is.null(fit$upper)) fit$upper <- fit$upper[keep]
	    if (!is.null(fit$lower)) fit$lower <- fit$lower[keep]
	    }
	}
    fit
    }

# SCCS  @(#)survobrien.s	4.3 04/14/92
#
# The test for survival proposed by Peter O'Brien
#
survobrien <- function(formula, data= sys.frame(sys.parent())) {
    m <- model.frame(formula, data, na.action= function(x) x )
    n <- nrow(m)
    Terms <- attr(m, 'terms')

    y <- model.extract(m, 'response')
    if (!inherits(y, "Surv")) stop ("Response must be a survival object")
    if (attr(y, 'type') != 'right') stop("Can only handle right censored data")

    # Figure out which are the continuous predictor variables
    m.name <- names(m)
    temp <- match(attr(Terms, 'term.labels'), m.name)
    cont <- NULL
    for (i in temp) {if (!is.factor(m[[i]])) cont <- c(cont, i)}
    if (is.null(cont)) stop ("No continuous variables to modify")

    keepers <- rep(T, length(m))  #The ones to be kept "as is"
    keepers[cont] <- F
    keepers[as.numeric(attr(Terms, 'response'))] <- F
    ord <- order(y[,1])
    x <- as.matrix(m[ord, cont])
    time <- y[ord,1]
    status <- y[ord,2]
    nvar <- length(cont)

    nline <- 0
    for (i in unique(time[status==1])) nline <- nline + sum(time >=i)
    start <- stop <- event <- double(nline)
    xx <- matrix(double(nline*nvar), ncol=nvar)
    ltime <- 0
    j<- 1
    keep.index <- NULL
    for (i in unique(time[status==1])) {
	who <- (time >=i)
	nrisk <- sum(who)
	if (nrisk<2) break
	temp <- apply(x[who,,drop=F], 2, rank)
	temp <- (2*temp -1)/ (2* nrisk)   #percentiles
	logit<- log(temp/(1-temp))           #logits
	deaths <- (status[who]==1 & time[who]==i)

	k <- seq(from=j, length=nrisk)
	start[k] <- ltime
	stop[k] <-  i
	event[k] <- deaths
	xx[k,] <- logit
	j <- j + nrisk
	ltime <- i
	keep.index <- c(keep.index, ord[who])
	}

    if (any(keepers)) {
	 temp <- list(m[keep.index, keepers], start, stop, event, xx)
	 names(temp) <- c(m.name[keepers], "start", "stop", "event",
				m.name[cont])
	 }
    else {
	temp <- list(start, stop, event, xx)
	names(temp) <- c(m.name[keepers], "start", "stop", "event",
				m.name[cont])
	}
    temp
    }
#SCCS @(#)survreg.control.s	4.1 11/19/92
survreg.control <- function(maxiter=30, rel.tolerance=1e-5, failure=1)
    list(maxiter = maxiter, rel.tolerance = rel.tolerance, failure =failure)
# SCCS @(#)survreg.distributions.s	4.3 11/19/92
#
# Create the survreg.distributions object
#
###<TSL> can't use quotes for names in list()
survreg.distributions <- list(
extreme = list(
    name = "Extreme value",
    parms = list(scale=1),
    logconcave = T,
    init  = function(x,fixed, init) {
			if (!is.null(fixed$scale))
				temp<- c(log(fixed$scale), 1)
			else if(!is.null(init$scale))
				temp<- c(log(init$scale), 0)
			else    temp<- c(log(mean(x^2)/1.3)/2 ,0)
			matrix(temp, nrow=1, dimnames=list("Log(scale)", NULL))
			},
    deviance= function(y, parms, loglik) {
			status <- y[,ncol(y)]
			scale <- exp(parms[1])
			if (any(status==3)) {
			    temp <- ifelse(status==3,(y[,2] - y[,1])/scale,1)
			    temp2 <- temp/(exp(temp)-1)
			    temp3 <- log(exp(-temp2) - exp(-temp2*exp(temp)))
			    best <- ifelse(status==1, -(1+log(scale)),
				    ifelse(status==3, temp3, 0))
			    }
			else best <- ifelse(status==1, -(1+log(scale)), 0)
			2*(best-loglik)
			},
    print = function(parms, fixed)
		if (fixed)
		     paste("Dispersion (scale) fixed at", format(exp(parms)))
		else paste("Dispersion (scale) =",
					       format(exp(parms)))
    ),

logistic = list(
    name  = "Logistic",
    parms = list(scale=1),
    logconcave = T,
    init  = function(x,fixed, init) {
			if (!is.null(fixed$scale))
				temp<- c(log(fixed$scale), 1)
			else if(!is.null(init$scale))
				temp<- c(log(init$scale), 0)
			else    temp<- c(log(mean(x^2)/2)/2 ,0)
			matrix(temp, nrow=1, dimnames=list("Log(scale)", NULL))
			},
    deviance= function(y, parms, loglik) {
			status <- y[,ncol(y)]
			scale <- exp(parms[1])
			if (any(status==3)) {
			    temp <- ifelse(status==3,(y[,2] - y[,1])/scale,1)
			    temp <- (y[,2] - y[,1])/scale
			    temp2 <- exp(temp/2)
			    temp3 <- log((temp2-1)/(temp2+1))
			    best <- ifelse(status==1, -log(4*scale),
				    ifelse(status==3, temp3, 0))
			    }
			else best <- ifelse(status==1, -log(4*scale), 0)
			2*(best-loglik)
			},
    print = function(parms, fixed)
		if (fixed)
		     paste("Dispersion (scale) fixed at", format(exp(parms)))
		else paste("Dispersion (scale) est =",
					       format(exp(parms)))
    ),

gaussian = list(
    name  = "Gaussian",
    parms = list(scale=1),
    logconcave = T,
    init  = function(x,fixed, init) {
			if (!is.null(fixed$scale))
				temp<- c(log(fixed$scale), 1)
			else if(!is.null(init$scale))
				temp<- c(log(init$scale), 0)
			else    temp<- c(log(mean(x^2))/2 ,0)
			matrix(temp, nrow=1, dimnames=list("Log(scale)", NULL))
			},
    deviance= function(y, parms, loglik) {
			status <- y[,ncol(y)]
			scale <- exp(parms[1])
			temp <-  ifelse(status==3, (y[,2] - y[,1])/scale, 1)
			temp2 <- 2*(pnorm(temp/2) -.5)
			best <- ifelse(status==1, -log(sqrt(2*pi)*scale),
				ifelse(status==3, temp2, 0))
			2*(best-loglik)
			},
    print = function(parms, fixed)
		if (fixed)
		     paste("Dispersion (scale) fixed at", format(exp(parms)))
		else paste("Dispersion (scale) =",
					       format(exp(parms)))
    ),

t = list(
    name  = "Student-t",
    parms = list(scale=1, df=2),
    logconcave = F,
    init  = function(x,fixed, init) {
			if (!is.null(fixed$scale))
				temp<- c(log(fixed$scale), 1)
			else if(!is.null(init$scale))
				temp<- c(log(init$scale), 0)
			else    temp<- c(log(mean(x^2))/2 ,0)

			if (!is.null(fixed$df))
				temp<- c(temp, fixed$df, 1)
			else if(!is.null(init$df))
				temp<- c(temp, init$df, 0)
			else    {
			    k <- mean(x^4)/ (3*temp[1])
			    df<- (3*k + sqrt(8*k + k^2))/(16*k*(k-1))
			    temp<- c(temp, df, 0)
			    }

			matrix(temp, nrow=2, byrow=T,
			    dimnames=list(c("Log(scale)", "df"),  NULL))
			},
    deviance= function(y, parms, loglik) {
			status <- y[,ncol(y)]
			scale <- exp(parms[1])
			df <- parms[2]
			temp <-  ifelse(status==3, (y[,2] - y[,1])/scale, 1)
			temp2 <- 2*(pt(temp/2, df) -.5)
			temp3 <- lgamma((df+1)/2) -
				    (lgamma(df/2) + .5*log(pi*df*scale^2))
			best <- ifelse(status==1, temp3,
				ifelse(status==3, temp2, 0))
			2*(best-loglik)
			},
    print = function(parms, fixed) {
	    tt <- if (fixed[1])
		     paste("Dispersion (scale) fixed at", format(exp(parms[1])))
		else paste("Dispersion (scale) =",
					       format(exp(parms[1])))
	    if (fixed[2]) paste(tt, ", df fixed at", format(parms[2]))
	    else tt
	    }
    )
 )
#SCCS @(#)survreg.fit.s	4.8 12/29/92
#
# This handles the one parameter distributions-- extreme, logistic,
#       gaussian, and cauchy.
# The parent routine survreg can't allow Cauchy.  It gives negative weights
#   when we try to fit it into the IRLS model, as Cauchy is not log-covex.
#
survreg.fit<- function(x, y, offset, init, controlvals, dist, fixed,
				nullmodel=T) {

    iter.max <- controlvals$maxiter
    eps <- controlvals$rel.tol

    if (!is.matrix(x)) stop("Invalid X matrix ")
    n <- nrow(x)
    nvar <- ncol(x)
    ny <- ncol(y)
    if (is.null(offset)) offset <- rep(0,n)

    sd <- survreg.distributions[[dist]]
    if (is.null(sd)) stop ("Unrecognized distribution")
    dnum <- match(dist, c("extreme", "logistic", "gaussian", "cauchy"))
    if (is.na(dnum)) stop ("Unknown distribution given")

    ncol.deriv <- if (any(y[,ny]==3)) 6  else 3
    nvar2 <- nvar + length(sd$parms) - length(fixed)
    yy <- ifelse(y[,ny]!=3, y[,1], (y[,1]+y[,2])/2 )

    if (is.numeric(init)) {
	if (length(init) != nvar2) stop("Wrong length for initial parameters")
	eta <- x %*% init[1:nvar]
	tfix <- sd$init(yy - eta, fixed, init)
	}
    else {
	if (is.null(eta <- init$eta))  eta <- mean(yy)
	else if (length(eta) != n) stop ("Wrong length for init$eta")

	# Do the 'glim' method of finding an initial value of coef,
	tfix <- sd$init(yy - (eta+offset), init, fixed)
	deriv <- .C("survreg_g",
		       as.integer(n),
		       y,
		       as.integer(ny),
		       as.double(eta + offset),
		       coef= as.double(c(0,tfix[,1])),
		       deriv = matrix(double(n * 3),nrow=n),
		       as.integer(3),
		       as.integer(dnum))$deriv
	wt <-  -1*deriv[,3]
	coef <- solve(t(x)%*% (wt*x), c((wt*eta + deriv[,2])%*% x))
	eta <- x %*% coef  + offset
	tfix <- sd$init(yy-eta, fixed, init)
	init <- c(coef, tfix[tfix[,2]==0,1])
	}

    fit <- .C("survreg",
		   iter = as.integer(iter.max),
		   as.integer(n),
		   as.integer(nvar),
		   y,
		   as.integer(ny),
		   x,
		   as.double(offset),
		   coef= as.double(init),
		   as.integer(nrow(tfix)),
		   tfix,
		   double(3*(nvar2) + 2*n),
		   imat= matrix(0, nvar2, nvar2),
		   loglik=double(2),
		   flag=integer(1),
		   as.double(eps),
		   deriv = matrix(double(ncol.deriv*n),nrow=n),
		   as.integer(dnum))

    if (iter.max >1 && fit$flag== -1) {
	if (controlvals$failure==1)
	       warning("Ran out of iterations and did not converge")
	else if (controlvals$failure==2)
	       return("Ran out of iterations and did not converge")
	}

    temp <- dimnames(x)[[2]]
    if (is.null(temp)) temp <- paste("x", 1:ncol(x))
    temp <- c(temp, (dimnames(tfix)[[1]])[tfix[,2]==0])
    dimnames(fit$imat) <- list(temp, temp)
    names(fit$coef) <- temp
    parms <- tfix[,1]
    parms[tfix[,2]==0] <- fit$coef[-(1:nvar)]

    if (!nullmodel)
	c(fit[c("iter", "coef", "imat", "loglik", "flag", "deriv")],
		list(parms=parms, fixed=tfix[,2]==1))
    else {
	init <- c(mean(x%*% fit$coef[1:nvar]), fit$coef[-(1:nvar)])
	temp <- cbind(parms, 1)     # "nail down" extras
	nfit <- .C("survreg",
		   iter = as.integer(iter.max),
		   as.integer(n),
		   nvar= as.integer(1),
		   y,
		   as.integer(ny),
		   x= rep(1,n),
		   as.double(offset),
		   coef= as.double(init),
		   as.integer(nrow(tfix)),
		   temp,
		   double(3*(nvar) + 2*n),
		   imat= matrix(0, nvar, nvar),
		   loglik=double(2),
		   flag=integer(1),
		   as.double(eps),
		   deriv = matrix(double(ncol.deriv*n),nrow=n),
		   as.integer(dnum))

	c(fit[c("iter", "coef", "imat", "loglik", "flag", "deriv")],
	       list(ndev=nfit$loglik, ncoef=nfit$coef, parms=parms,
		    fixed=tfix[,2]==1))
	}
    }
#SCCS @(#)survreg.s	4.15 02/04/95
survreg <- function(formula=formula(data), data=sys.frame(sys.parent()),
	subset, na.action,
	link='log',
	dist=c("extreme", "logistic", "gaussian", "exponential",
	       "rayleigh"),
	init=NULL,  fixed=list(), control,
	model=F, x=F, y=T, ...) {

    call <- match.call()
    m <- match.call(expand=F)
    m$dist <- m$link <- m$model <- m$x <- m$y <- m$... <-  NULL
    m$start <- m$fixed <- m$control <- m$init <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.frame(sys.parent()))
    Terms <- attr(m, 'terms')

    dist <- match.arg(dist)
##    lnames <- dimnames(glm.links)[[2]]
##    link <- pmatch(link, lnames, 0)
    link<-make.link(link)
##    if (link==0) stop("Invalid link function")
##    else link <- lnames[link]

    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    X <- model.matrix(Terms, m)
    n <- nrow(X)
    nvar <- ncol(X)
    offset<- attr(Terms, "offset")
    if (!is.null(offset)) offset <- as.numeric(m[[offset]])
    else                  offset <- rep(0, n)

    type <- attr(Y, "type")
##    linkfun <- glm.links["link", link][[1]]
    linkfun<-link$linkfun
    if (type== 'counting') stop ("Invalid survival type")
    else if (type=='interval') {
	if (any(Y[,3]==3))
	     Y <- cbind(linkfun(Y[,1:2]), Y[,3])
	else Y <- cbind(linkfun(Y[,1]), Y[,3])
	}
    else if (type=='left')
	     Y <- cbind(linkfun(Y[,1]), 2-Y[,2])
    else     Y <- cbind(linkfun(Y[,1]), Y[,2])

    controlvals <- survreg.control()
    if (!missing(control)) 
	controlvals[names(control)] <- control

    if( dist=="exponential") {
	fixed$scale <- 1
	dist <- 'extreme'
	}
    else if (dist=="rayleigh") {
	fixed$scale <- .5
	dist <- 'extreme'
	}

    sd <- survreg.distributions[[dist]]
    if (length(fixed)>0) {
	ifix <- match(names(fixed), names(sd$parms), nomatch=0)
	if (any(ifix==0))
	    stop (paste("Parameter(s)", paste(names(fixed)[ifix==0]),
			"in the fixed list not valid for this dist"))
	}
    if (is.list(init) && length(init)>0) {
	ifix <- match(names(init), c('eta',names(sd$parms)), nomatch=0)
	if (any(ifix==0))
	    stop (paste("Parameter(s)", paste(names(init)[ifix==0]),
			"in the init list not valid for this dist"))
	}


    sfit <- survreg.fit(X, Y, offset, init=init, controlvals=controlvals,
			dist= dist, fixed=fixed)
    if (is.character(sfit))  fit <- list(fail=sfit)  #error message
    else {
	# There may be more clever ways to do this, but ....
	#  In order to make it look like IRLS, and get all the components
	#  that I need for glm inheritance, do one step of weighted least
	#  squares.
	eta <- c(X %*% sfit$coef[1:nvar]) + offset
	wt<- -sfit$deriv[,3]
###	fit <- lm.wfit(X, eta + sfit$deriv[,2]/wt, wt, "qr", ...)
	fit <- lm.wfit(X, eta + sfit$deriv[,2]/wt, wt)
###	ifun <- glm.links['inverse',link][[1]]
	ifun <-link$linkinv
	fit$fitted.values <- ifun(fit$fitted.values)
	fit$family <- list(name=dist, link=link, "")
	fit$linear.predictors <- eta
	fit$iter <- sfit$iter
	fit$parms <- sfit$parms
	fit$df.residual <- fit$df.residual - sum(!sfit$fixed)

	# If singular.ok=T, there may be NA coefs.  The var matrix should
	#   be an inversion of the "non NA" portion.
	var <- 0*sfit$imat
	good <- c(!is.na(fit$coef), rep(T, ncol(var)-nvar))
	var[good,good] <- solve(qr(sfit$imat[good,good], tol=1e-12))
	fit$var <- var
	fit$fixed <- sfit$fixed
	dev <- sd$deviance(Y, fit$parms, sfit$deriv[,1])
	fit$dresiduals <- sign(fit$residuals)*sqrt(dev)
	fit$deviance <- sum(dev)
	fit$null.deviance <- fit$deviance +2*(sfit$loglik[2]- sfit$ndev[2])
	fit$loglik <- c(sfit$ndev[2], sfit$loglik[2])
	}

    na.action <- attr(m, "na.action")
    if (length(na.action)) fit$na.action <- na.action
    attr(fit, "class") <-  c("survreg", "glm", "lm")
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- call
    if (model) fit$model <- m
    if (x)     fit$x <- X
    if (y)     fit$y <- Y
    fit
    }
# @(#)survsum.s	1.2 09/25/92
# written by Mark Dietrich
survsum <- function	
		(formula,data=sys.frame(sys.parent()),sptms=NULL,xlim,tlines=T,log=F,
		xscale=1,yscale=100,mark.time=F,mark=3,cex=1,xlab="Time",
		ylab="Survival (%)",lgd="bl",ttl="K-M Survival",...) {
#
	ltms <- length (sptms)			##number of specified times
#
#------------------------------------------------------------------------------
#
	if (ltms >4){
        	stop("Maximum number of specified times is four.")}
	if (ltms > 0){						##warnings
		if( any (sptms < 0))
			stop ("specified times must be positive")}
#
#------------------------------------------------------------------------------
#################
#total statistics#
##################
#
	fit <- survfit (formula,data)			##survfit object
#
	n.pts <- summary (fit,times=0,print.it=F)$n.risk
					      ##number of points in each group
	strat <- fit$strata		      ##
	gp.name <- names(strat)	      ## the group names
	n <- length (n.pts)		      ##the number of groups
#
	if (n > 6) {					##too many groups
		stop("Maximum number of groups is 6")}
#
	code <- fit$n.event				##coded events by time
	stemp <- rep (1:n,strat)
		events <- 1:n		     ##number of events in each group 
        		for (i in (1:n))
                		events[i] <- sum (code [stemp ==i])
#
#------------------------------------------------------------------------------
###############
#survival plot#
###############
#
	par (fig = c(0,1,.10,.75))
	par (err=-1)			       ##supress out-of-bounds msgs
#
	frame ()
	if (missing(xlim)){		##conditional: no xlim specified
#
		if (log) {			##conditional: log=True
#
	ends <- plot (fit,lty=c(1,2,7,8,3,5),xaxs="i",xscale=xscale,log=T,
	,xlab=xlab,ylab=ylab,mark.time=mark.time,mark=mark,cex=cex,las=1,...)
								         ##Plot
		} else {			##conditional: log=False
#
	ends <- plot (fit,lty=c(1,2,7,8,3,5),xaxs="i",yaxs="i",
	xscale=xscale,yscale=yscale,xlab=xlab,ylab=ylab,mark.time=mark.time,
	mark=mark,cex=cex,ylim=c(0,(yscale+(yscale/10))),las=1,...)	##Plot
#
		}
	xlim <- c(0,max(ends$x))       	##supply xlim value needed below
#
	} else {			##conditional: xlim is specified
#
		if (log) {			##conditional: log=True
#
	ends <- plot (fit,lty=c(1,2,7,8,3,5),xaxs="i",xlim=xlim,
	log=T,xscale=xscale,xlab=xlab,ylab=ylab,mark.time=mark.time,mark=mark,
	cex=cex,las=1,...)     						##Plot
#
		} else {			##conditional: log=False
#
	ends <- plot (fit,lty=c(1,2,7,8,3,5),xaxs="i",yaxs="i",
	xlim=xlim,xlab=xlab,ylab=ylab,mark.time=mark.time,mark=mark,cex=cex,
	xscale=xscale,yscale=yscale,ylim=c(0,(yscale+(yscale/10))),las=1,...)
									##Plot
			}}
#
#
	small <- xlim[1]			##minimum x value on plot 
	big <- xlim[2]				##maximum x value on plot
#	
	par (err=0)				##error msgs resumed
#
#------------------------------------------------------------------------------
#################
#specified times#
#################
#
	vec <- rep(NA,4)		##vector of NA's:used if ltms=0
	vec[1:ltms] <- sptms		##NA's replaced with specified times
	t1 <- vec[1]			##variables assigned timepoints
	t2 <- vec[2]
	t3 <- vec[3]
	t4 <- vec[4]
#
#-----------------------------------------------------------------------------
################
#vertical lines#
################
#
        if (tlines){				##conditional: tlines=True
		lin.ok <- vec [!is.na(vec)]	##times that are not NA
#
	if (ltms != 0){
	if (any (lin.ok < small) | any (lin.ok > big))
						##conditional:	times <
						##             xmin or > xmax
		stop ("a specified (sptms) time point falls outside the plot region.")
#
		abline (v=c(lin.ok),lty=4)		 ##vertical lines
                axis (3,at=c(lin.ok),labels=c(lin.ok))}} ##axis labels at times
#
#------------------------------------------------------------------------------
########
#legend#
########
#
	par (fig=c(0,1,0,1))
	a <- 1:100
	b <- 1:100
#
	plot (a,b,type="n",xlab="",ylab="",bty="n",axes=F)
#
	if (lgd=="bl") {		##conditional: legend at bottom-right
#
	if (n == 6) {
		legend (0,30,legend=names (fit$strata),
							lty=c(1,2,7,8,3,5))
			} else {
	if (n == 5) {
		legend (0,28,legend=names (fit$strata),
							lty=c(1,2,7,8,3,5))
			} else {
	if (n == 4) {
		legend (0,26,legend=names (fit$strata),
							lty=c(1,2,7,8,3,5))
			} else {
	if (n == 3) {
		legend (0,24,legend=names (fit$strata),
							lty=c(1,2,7,8,3,5))
			} else {
	if (n == 2) {
		legend (0,22,legend=names (fit$strata),
							lty=c(1,2,7,8,3,5))
			} else {
		legend (0,20,legend=names (fit$strata),
							lty=c(1,2,7,8,3,5))
}}}}}
        } else {
#
	if (lgd =="tr") {		##conditional: legend at top-right
#
		legend (75,69,legend=names (fit$strata),
		lty=c(1,2,7,8,3,5))
#
        } else {
#
	if (lgd == "n") 			##conditional: no legend 
		{	}}}
#
	par (fig = c(0,1,0,1))
#
	if (lgd == "under"){		##conditional: legend under plot
		legend (75,4,legend=names (fit$strata),
		lty=c(1,2,7,8,3,5))}
#
#------------------------------------------------------------------------------
#####################
#test for difference#
#####################
#
	sdif <- survdiff(eval(fit$call$formula), eval (fit$call$data))
					##survdiff function
        chisq <- round(sdif$chisq,2)	##chisquare statistic
        df <- length (sdif$n) - 1	##degrees of freedom
        p <- round(1 - pchisq(sdif$chisq, df),4)	##p-value
                if (p < .0001) (p <- ".0001")		
#
	text (0,-5,paste("LOGRANK TEST (all curves equal): Chi-Square = "
	,chisq,"   ","df = ",df,"   ","p = ",p),adj=0)
#
	mtext(date(),side=1,outer=T,adj=1)		##date
#
#------------------------------------------------------------------------------
#############################
#printing on graphics device#
#############################
#
	ysfull <- c(80,72,64,56,48,40)          ##y-values
        ys <- ysfull[1:n]
#
	par (fig = c(0,1,.5,1))
#
# jds -- switch ttl and "K-M Survival"
	if (!missing(ttl)) 
		mtext("K-M Survival",side=3,outer=T,adj=1)
	plot (a,b,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",main=ttl)
#
	x1 <- c(rep(-5,n+1),-5,20)
	y1 <- c(90,ys,90,105)
	labs1 <- c("Group",gp.name,"_________________________________________",
			"Totals")
	text (x1,y1,labs1,adj=0)
#
	x2 <- c(20,rep(18,n),32,rep(30,n))
	y2 <- c(90,ys,90,ys)
	labs2 <- c("No.Pts.",n.pts,"No.Events",events)
	text (x2,y2,labs2,adj=1)
#
#
##########################
#specific time statistics#
##########################
#
	gt1 <- gt2 <- gt3 <- gt4 <- NA		##return value dummy variables

	if (!is.na (t1)) {			##conditional: t1 is not NA
#
	text (38,105,"Estimated survival % (SE,n) at time (t)",adj=0)
#
#################
#define function#
#################
#
	group <- function (ti,x,current.fit,m,gpn,scale,endsx) {
#
        	        if (ti > max(endsx)){	##conditional: time > xmax
#
	        text (x,90,c(paste("t=", ti),"___________"),adj=0)
                text(x,80,"no data",adj=0)
        	        } else {		##conditional: time < xmax
#
	        mat <- summary.survfit (current.fit,times=ti,print.it=F,
					scale=scale)
#	
		gps <- mat$strata	##group names used at time

	        bigm <- matrix(rep(NA,m*3),ncol=3)
#
	        dimnames (bigm) <- list (gpn,c("percs","se","n"))
#
        	percs <- format (round (mat$surv*100,1))  ##survival percentage
        	sters <- format(round (mat$std.err*100,1))     ##standard error
        	nrisk <- format(mat$n.risk)		       ##no. at risk
#
		bigm [as.character(gps),] <- c(percs,sters,nrisk)

        	ysfull <- c(80,72,64,56,48,40)          ##y-values
        	ys <- ysfull[1:m]
#
        	text (x,90,c(paste("t=", ti),"___________"),adj=0)
        	text (rep(x,m),ys,paste(bigm[1:m,1],"(",bigm[1:m,2],",",
		bigm[1:m,3],")"),adj=0)
		list (bigm = bigm)
		}}
#
		f1 <- group (t1,35,fit,n,gp.name,xscale,ends$x)
		gt1 <- f1$bigm
		gt2 <- gt3 <- gt4 <- NA
#
	if (!is.na (t2)) {			##conditional: t2 is not NA
		f2 <- group (t2,52,fit,n,gp.name,xscale,ends$x)
		gt2 <- f2$bigm
		gt3 <- gt4 <- NA
#
	if (!is.na (t3)) {			##conditional: t3 is not NA
		f3 <- group (t3,69,fit,n,gp.name,xscale,ends$x)
		gt3 <- f3$bigm
		gt4 <- NA
#
	if (!is.na (t4)) {			##conditional: t4 is not NA
		f4 <- group (t4,86,fit,n,gp.name,xscale,ends$x)
		gt4 <- f4$bigm
}}}}
invisible(list(no.pts=n.pts,no.events=events,chisq=chisq,p=p,t1=gt1,
		t2=gt2,t3=gt3,t4=gt4 ))			##return values
}
















#SCCS @(#)tcut.s	4.2 10/16/94
tcut <-  function (x, breaks, labels, scale=1){
    if(length(breaks) == 1) {
	if(breaks < 1)
		stop("Must specify at least one interval")
	if(missing(labels))
		labels <- paste("Range", seq(length = breaks))
	else if(length(labels) != breaks)
		stop("Number of labels must equal number of intervals")
	r <- range(x[!is.na(x)])
	r[is.na(r)] <- 1
	if((d <- diff(r)) == 0) {
		r[2] <- r[1] + 1
		d <- 1
	    }
	breaks <- seq(r[1] - 0.01 * d, r[2] + 0.01 * d, length = breaks +1)
	}
    else {
	if(is.na(adb <- all(diff(breaks) >= 0)) || !adb)
	   stop("breaks must be given in ascending order and contain no NA's")
	if(missing(labels))
	    labels <- paste(format(breaks[ - length(breaks)]),
			"+ thru ", format(breaks[-1]), sep = "")
	else if(length(labels) != length(breaks) - 1)
	   stop("Number of labels must be 1 less than number of break points")
	}

    structure(x*scale, cutpoints=breaks*scale, labels=labels, class='tcut')
    }

"[.tcut" <- function(x,i) {
    atts <- attributes(x)
    class(x) <- NULL
    x <- x[i]
    attributes(x) <- atts
    x
    }
#SCCS @(#)untangle.specials.s	4.1 07/10/92
untangle.specials <- function(tt, special, order=1) {
    #
    # There was a change in the specials, so action depends on your release
    #   of S
    #
    spc <- attr(tt, 'specials')[[special]]
    if (length(spc)==0)
	return(vars=character(0), terms=numeric(0))

###<TSL>    facs <- attr(tt, 'factor')
    facs <- attr(tt, 'factors')
###</TSL>
    fname <- dimnames(facs)

    if ((attr(terms(y~ zed(x), specials='zed'), 'specials'))$zed ==1) {
	# old style
	if (any(order>1))
	   warning("Can't select specials based on the order of the terms")
	list(vars=(fname[[2]])[spc],  terms=spc)
	}
    else {
	ff <- apply(facs[spc,,drop=F], 2, sum)
	list(vars= (fname[[1]])[spc],
	     terms= seq(ff)[ff & match(attr(tt, 'order'), order, nomatch=0)])
	}
    }
#SCCS 04/14/92 @(#)model.newframe.s	4.3
# This function is called if you want to get a new data frame,
#   usually for prediction.  It's main problem is to "glue" any
#   transform specific information back onto the formula, so that
#   data dependent transforms work as they used to.
# It only works if the data dependent functions are not inside another one,
#   so  sqrt(age - min(age)) is out of luck.  It also only works for those
#   transforms that support it by adding data dependent info as an attribute
#   of their output.
# If you know this isn't so, then safe=T uses a method that is much longer,
#   but is guarranteed to work, see predict.gam

model.newframe <- function(object, newdata, safe=F, response=F, ...) {
    if (inherits(object, 'terms'))  Terms <- object
    else {
	Terms <- object$terms
	if (!inherits(Terms, 'terms'))
	    stop ("Invalid terms component of object")
	}
    offset <- attr(Terms, 'offset')

    # First, is newdata just a list of numbers?
    if (is.numeric(newdata)) {
##	nvar <- length(Terms) + length(offset)
        nvar <- length(attr(Terms,"term.labels")) +length(offset)
	if (length(newdata)>1  || newdata!=floor(newdata)  || newdata<0){ #It's not just a frame number
	    if (is.matrix(newdata) && ncol(newdata) == nvar)
		   return(newdata)
	    else if (length(newdata) == nvar)
		   return(matrix(newdata,1,nvar))
	    else stop("Argument \"newdata\" cannot be coerced to an appropriate model matrix")
	    }
	}

    # newdata is a list, data frame, or frame number
    if (!safe) {
	#augment the arguments with extra parameters
	  #someday
	if (!response) Terms <- delete.response(Terms)
	model.frame(Terms, newdata, ...)
	}
    else {
	#Do a safe call, by building up a brand new model frame
	Call <- object$call
	Call[[1]] <- as.name("model.frame")
	Call$formula <- terms.inner(formula(object))
   #might need to tack on the response here!
	if (response) stop("Not implimented yet for safe=T, response=T")
	Call$na.action <- function(x)  x
	Call <- Call[match(c("", "formula", "data", "subset", "na.action"),
	    names(Call), 0)]
	data <- eval(Call)
	attr(data, "terms") <- NULL
	Call$subset <- NULL
	Call$data <- substitute(newdata)
	newdata <- eval(Call)
	attr(newdata, "terms") <- NULL
	d2 <- dim(newdata)
	if(d2[1] < 1)
	    stop("0 rows in newdata")
	d1 <- dim(data)
	if(d1[2] != d2[2])  #newdata missing some variables
	    data <- data[, names(newdata), drop = F]
	data[seq(d1[1] + 1, d1[1] + d2[1]),  ] <- newdata  #rbind the new on
	attr(data, "row.names") <- c(rep("OLD DATA",d1[1]), row.names(newdata))
	#Now compute the combined model frame, excluding the response
	na.action <- eval(object$call$na.action)
	Terms <- object$terms
	Terms <- delete.response(Terms)
	model.frame(Terms, data, na.action = na.action)
	}
    }
