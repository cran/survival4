.First.lib <- function(lib, pkg) {
	library.dynam(pkg, lib.loc=lib)
}	

if (exists("autoload")){
	autoload("survexp.us","ratetables")
	autoload("survexp.uswhite","ratetables")
	autoload("survexp.usr","ratetables")
	autoload("survexp.az","ratetables")
	autoload("survexp.azr","ratetables")
	autoload("survexp.fl","ratetables")
	autoload("survexp.flr","ratetables")
	autoload("survexp.mn","ratetables")
	autoload("survexp.mnwhite","ratetables")
	autoload("survexp.wnc","ratetables")
	}




"MyUseMethod" <-
function (generic, classobj = NULL) 
{ ## This probably isn't needed any more but I haven't got around to
  ## checking (TSL 25/8/97)
        call <- sys.call(sys.parent())
	if (is.null(classobj)) 
		classobj <- eval(call[[2]], sys.frame(sys.parent()))
	classlist <- class(classobj)
	classlist <- c(classlist, "default")
	while (!exists(paste(generic, classlist[[1]], sep = "."
	),mode="function",inherits=TRUE) && length(classlist) > 1) classlist <- classlist[-1]
        methodname<-paste(generic, classlist[[1]], sep = "."
	)
	if (!exists(methodname,mode="function",inherits=TRUE))
		stop("No method found")
	call[[1]] <- as.name(methodname)
	eval(call, sys.frame(-2))
}

all.equal <- function(x,y){
	 as.logical(prod(as.numeric((abs(x-y)<.Machine$double.eps) || (is.na(x) && is.na(y)))))
}


 sort.list<-function(...) order(...)

"[.terms" <-
function (termobj, i) 
{
        resp <- if (attr(termobj, "response")) 
                termobj[[2]]
        else NULL
        newformula <- attr(termobj, "term.labels")[i]
        if (length(newformula) == 0) 
                newformula <- 1
        newformula <- reformulate(newformula, resp)
        terms(newformula, specials = names(attr(termobj, "specials")))
}

is.category <- function(x) inherits(x,"factor") || is.factor(x)

surv4.prmatrix<-function(mat,rowlab=NULL,collab=NULL,...){
  if (is.null(rowlab)) rowlab<-dimnames(mat)[[1]]
  if (is.null(collab)) collab<-dimnames(mat)[[2]]
  prmatrix(mat,rowlab=rowlab,collab=collab,...)
}
