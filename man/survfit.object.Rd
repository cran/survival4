\name{survfit.object}
\title{
Parametric Survival Model Object
}
\description{
This class of objects is returned by the `survfit' class of functions
to represent a fitted survival curve.


Objects of this class have methods for the functions `print',
`summary', `plot', `points' and `lines'.
}
\section{COMPONENTS}{
The following components must be included in a legitimate \code{survfit} object.
}
\arguments{
\item{time}{
the time points at which the curve has a step.
}
\item{n.risk}{
the number of subjects at risk at t.
}
\item{n.event}{
the number of events that occur at time t.
}
\item{surv}{
the estimate of survival at time t+0.
This may be a vector or a matrix.
}
\item{strata}{
if there are multiple curves, this component gives the number of elements of
of the \code{time} etc. vectors corresponding to the first curve, the second curve,
and so on.  The names of the elements are labels for the curves.
}
\item{std.err}{
the standard error of the cumulative hazard or -log(survival).
}
\item{uppper}{
upper confidence limit for the survival curve.
}
\item{lower}{
lower confidence limit for the survival curve.
}
\item{conf.type}{
the approximation used to compute the confidence limits.
}
\item{conf.int}{
the level of the confident limits, e.g. 90 or 95\%.
}
\item{na.action}{
the returned value from the na.action function, if any.  It will be used
in the printout of the curve, e.g., the number of observations deleted due
to missing values.
}
\item{call}{
an image of the call that produced the object.
}}
\seealso{
\code{\link{survfit}}, \code{\link{plot.survfit}}, \code{\link{summary.survfit}}.
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
