# " SCCS @(#)coxph.d	1.6 10/2/94)
\name{coxph}
\title{
Proportional Hazards Regression
}
\usage{
coxph(formula=formula(data), data=sys.parent(), subset, 
       na.action, weights, eps=0.0001, init,
       iter.max=10, method=c("efron","breslow","exact"),
       singular.ok=T, robust,
       model=F, x=F, y=T)
}
\arguments{
\item{formula}{
a formula object, with the response on the left of a ~ operator, and
the terms on the right.  The response must be a survival object as
returned by the Surv() function.
}
\item{data}{
a data.frame in which to interpret the variables named in
the formula, or in the subset and the weights argument.
}
\item{subset}{
expression saying that only a subset of the rows of the data
should be used in the fit.
}
\item{na.action}{
a missing-data filter function, applied to the model.frame, after any
subset argument has been used.  Default is options()$na.action.
}
\item{weights}{
case weights.
}
\item{eps}{
convergence criteria.  Iteration will continue until relative change
in log-likelihood is less than eps.  Default is .0001.
}
\item{init}{
vector of initial values of the iteration.  Default initial
value is zero for all variables.
}
\item{iter.max}{
maximum number of iterations to perform.  Default is 10.
}
\item{method}{
method for tie handling.  If there are no tied death times all the methods
are equivalent.
Nearly all Cox regression programs use the Breslow method by default, but
not this one.
The Efron approximation is used as the default here, as it is much more
accurate when dealing with tied death times, and is as efficient
computaionally.
The exact method computes the exact partial likelihood, which is equivalent
to a conditional logistic model.  If there are a large number of ties the
computational time will be excessive.
}
\item{singular.ok}{
If TRUE, the program will automatically skip over columns of the X matrix
that are linear combinations of earlier columns.  In this case the
coefficients for such columns will be NA, and the variance matrix will contain
zeros.  For ancillary calculations, such as the linear predictor, the missing
coefficients are treated as zeros.
}
\item{robust}{
if TRUE a robust variance estimate is returned.  Default is TRUE if the
model includes a cluster() operative, FALSE otherwise.
}
\item{model,x,y}{
flags to control what is returned.  If these are true, then the model
frame, the model matrix, and/or the response is returned as components
of the fitted model, with the same names as the flag arguments. 
}}
\value{
an object of class "coxph"
}
\section{Side Effects}{
Depending on the call, the predict, residuals, and survfit routines may
need to reconstruct the x matrix created by coxph.  Differences in the
environment, such as which data frames are attached or the value of
options()$contrasts, may cause this computation to fail or worse, to be
incorrect.  See the survival overview document for details.
}
\description{
The proportional hazards model is usually expressed in terms of a
single survival time value for each person, with possible censoring.
Anderson and Gill reformulated the same problem as a counting process;
as time marches onward we observe the events for a subject, rather
like watching a Geiger counter.
The data for a subject is presented as multiple rows or "observations", each
of which applies to an interval of observation (start, stop].
}
\section{CONVERGENCE}{
In certain data cases the actual MLE estimate of a
coefficient is infinity, e.g., a dichotomous variable where one of the
groups has no events.  When this happens the associated coefficient
grows at a steady pace and a race condition will exist in the fitting
routine: either the log likelihood converges, the information matrix
becomes effectively singular, an argument to exp becomes too large for
the computer hardware, or the maximum number of interactions is exceeded.
The routine attempts to detect when this has happened, not always
successfully.
}
\references{
Terry Therneau, author of local function.


P. Andersen and R. Gill. "Cox's regression model for
counting processes, a large sample study", Annals of Statistics, 
10:1100-1120, 1982.  


T.Therneau, P. Grambsch, and T.Fleming. "Martingale based residuals
for survival models", Biometrika, March 1990.
}
\seealso{
\code{\link{survfit}}, \code{\link{Surv}}, \code{\link{strata}}.
}
\examples{
# Create the simplest test data set
#
 test1 <- list(time=  c(4, 3,1,1,2,2,3),
                status=c(1,NA,1,0,1,1,0),
                x=     c(0, 2,1,1,1,0,0),
                sex=   c(0, 0,0,0,1,1,1))
 coxph( Surv(time, status) ~ x + strata(sex), test1)  #stratified model


#
# Create a simple data set for a time-dependent model
#
test2 <- list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
                stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
                event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
                x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) )


summary( coxph( Surv(start, stop, event) ~ x, test2))


}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
