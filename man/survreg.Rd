\name{survreg}
\title{
Regression for a parametric survival model
}
\usage{
survreg(formula, data=sys.parent(), subset, na.action,
link=c("log", "identity"),
dist=c("extreme", "logistic", "gaussian", "exponential"),
fixed, eps=0.0001, init, iter.max=10, model=F, x=F, y=F, \dots)
}
\arguments{
\item{formula}{
a formula expression as for other regression models.
See the documentation for \code{lm} and \code{formula} for details.
}
\item{data}{
optional data frame in which to interpret the variables occuring in the
formula.
}
\item{subset}{
subset of the observations to be used in the fit.
}
\item{na.action}{
function to be used to handle any NAs in the data.
}
\item{link}{
transformation to be used on the y variable.
}
\item{dist}{
assumed distribution for the transformed y variable.
}
\item{fixed}{
a list of fixed parameters, most often just the scale.
(When I implement the t-dist, it will include the degrees of freedom).
}
\item{eps}{
convergence criteria for the computation.  Iteration continues until the
relative change in log likelihood is less than eps.
}
\item{init}{
optional vector of initial values for the paramters.
}
\item{iter.max}{
maximum number of iterations to be performed.
}
\item{model}{
if TRUE, the model frame is returned.
}
\item{x}{
if TRUE, then the X matrix is returned.
}
\item{y}{
if TRUE, then the y vector (or survival times) is returned.
}
\item{...}{
all the optional arguments to lm, including \code{singular.ok}.
}}
\value{
an object of class \code{survreg} is returned, which inherits from class \code{glm}.
}
\section{Computation}{
  This routine is not as robust against nearly singular X matrices as lm();
the problem occurs when we explicitly invert the covariance matrix with
solve().  This can sometimes be solved by subtracting the mean from all
continuous covariates.
}
\examples{
data(ovarian)
survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='extreme',
		link='log', fixed=list(scale=1))   #Fit an exponential
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
