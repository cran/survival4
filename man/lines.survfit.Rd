\name{lines.survfit}
\title{
Add lines to a survival plot
}
\usage{
lines.survfit(x, type="s", mark=3, col=1, lty=1, lwd=1, mark.time=T, \dots)
}
\arguments{
\item{x}{
a survival object, generated either from the survfit or survexp functions.
}
\item{type}{
the line type, as described in \code{lines}.  The default is a step function
for survfit objects, and a connected line for survexp objects.
}
\item{mark,}{
vectors giving the mark symbol, color, line type and line width for the
added curves.
}
\item{mark.time}{
controls the labeling of the curves.  If False, no labeling is done.  If
True, then curves are marked at each censoing time.  If mark.time is a numeric
vector, then curves are marked at the specified time points.
}}
\section{Side Effects}{
one or more curves are added to the current plot.
}
\section{NOTE}{
Does not yet handle confidence intervals.
}
\seealso{
\code{\link{plot.survfit}}, \code{\link{survfit}}, \code{\link{survexp}}
}
\examples{
\dontrun{
data(jasa)
library(ratetables)
fit <- survfit(Surv(futime, fustat) ~ surgery, jasa)
ptime <- ifelse(jasa$fustat==0, jasa$futime, mdy.date(4,1,74)-jasa$accept.dt)
age <- jasa$accept.dt - jasa$birth.dt
efit <- survexp(ptime~ratetable(age=age, year=accept.dt, sex=1) + surgery,
		data=jasa, ratetable=survexp.us, conditional=F,times=0:5*300)
plot(fit, col=1:2, lty=1)
lines(efit, col=1:2, lty=2, mark='E', mark.time=1000, cex=1.5)
}
# Note: the 2 expected curves are overlap
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
