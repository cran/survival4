\name{survreg.object}
\title{
Parametric Survival Model Object
}
\description{
This class of objects is returned by the `survreg' function
to represent a fitted parametric survival model.
Class `survreg' inherits from class `glm', since it is fit by iterative
reweighted least squares; the object returned has all the components of a
weighted least squares object.


Objects of this class have methods for the functions `print',
`summary', `predict', and 'residuals'.
}
\section{COMPONENTS}{
The following components must be included in a legitimate \code{survreg} object.
The residuals, fitted values, coefficients and effects should be extracted
by the generic functions of the same name, rather than
by the \code{"$"} operator. 
}
\arguments{
\item{coefficients}{
the coefficients of the \code{linear.predictors}, which multiply  the
columns of the model
matrix.
It does not include the estimate of error (sigma).
The names of the coefficients are the names of the
single-degree-of-freedom effects (the columns of the
model matrix).
If the model is over-determined there will
be missing values in the coefficients corresponding to inestimable
coefficients.
}
\item{parms}{
the parameters of the model that are not coefficients of the X matrix.
The first of these will always be \code{log(sigma)}.
}
\item{fixed}{
a vector of the same length as \code{parms}, where 1 indicates a parameter that
was fixed at its starting value and was not part of the iteration.
}
\item{deviance}{
minus twice the difference between the maximized log-likelihood under the
fitted model and a saturated model.
Similar to the residual sum of squares.
}
\item{loglik}{
the log-likelihood for the final model.
}
\item{null.deviance}{
the deviance corresponding to the model with only an itercept term, and
with \code{parms} fixed at their final values.
}
\item{dresiduals}{
the deviance residuals.
}
\item{var}{
the final variance matrix, including both coefficients and free parameters.
}
\item{family}{
a 2 element character vector giving the name of the family and
the link; mainly for printing purposes.


The object will also have the components of an \code{glm} object:
\code{linear predictors}, \code{fitted.values}, \code{residuals},
\code{effects}, \code{R}, \code{rank}, \code{assign}, \code{contrasts}, \code{weights}, \code{iter},
\code{residuals}, \code{fitted.values}, \code{call}, \code{terms} and \code{formula}.
See \code{glm.object}.
}}
\seealso{
\code{\link{survreg}}, \code{\link{glm.object}}, \code{\link{lm.object}}.
}
\keyword{regression}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
