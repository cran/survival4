\name{summary.survfit}
\title{
Print a Survival Curve
}
\usage{
summary.survfit(fit, times, censored=F, scale=1)
}
\arguments{
\item{fit}{
output from a call to \code{survfit}.
}
\item{times}{
vector of the times;
the returned matrix will contain 1 row for each time.
This must be in increasing order and missing values are not allowed.
If censored=T, the default time vector contains all the unique times in
\code{fit.list},
otherwise the default time vector uses only the event (death) times.
}
\item{censored}{
logical flag: should the censoring times be included in the output?
This is ignored if the \code{times} argument is present.
}
\item{scale}{
rescale the survival time, e.g., if the input data to survfit were in
days, "scale=365" would scale the printout to years.
}
\value{
a list with the following components
}
\item{time}{
the timepoint on the curve.
}
\item{surv}{
the value of the survival curve at time t+0.
}
\item{n.risk}{
the number of subjects at risk at time t-0
(but see the comments on weights in the \code{survfit} help file).
}
\item{n.event}{
if the \code{times} argument is missing, then this column is the number of
events that occurred at time t.
Otherwise, it is the cumulative number of events that have occurred
since the last time listed until time t+0.
}
\item{std.dev}{
the standard deviation of the survival value.
}
\item{lower}{
lower confidence limits for the curve.
}
\item{upper}{
upper confidence limits for the curve.
}
\item{strata}{
the number of points per curve.  If there are multiple curves in this
result, then the \code{time}, \code{surv}, etc vectors will contain multiple
curves, pasted end to end, and length(surv)=sum(strata).  The names
of the strata vector are the labels for the curves.
}
\item{call}{
the statement used to create the \code{fit} object.
}
\item{na.action}{
passed through from \code{fit}, if present.
}}
\seealso{
\code{\link{survfit}}, '\code{\link{print.summary.survfit}}.
}
\examples{
data(ovarian)
summary( survfit( futime, fustat))
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
