\name{plot.cox.zph}
\title{
Graphical test of proportional hazards
}
\usage{
plot.cox.zph(x, resid=T, se=T, df=4, nsmo=40, var)
}
\arguments{
\item{x}{
result of the \code{cox.zph} function.
}
\item{resid}{
include the residuals on the plot, as well as the smooth fit.
}
\item{se}{
if true, confidence bands at 2 standard errors will be added.
}
\item{df}{
the degrees of freedom for the fitted natural spline.  A df value of
2 leads to a linear fit.
}
\item{nsmo}{
number of points used to plot the fitted spline.
}
\item{var}{
the set of variables for which plots are desired.  By default, plots are
produced in turn for each variable of a model.  Selection of a single variable
allows other features to be added to the plot, e.g., a horizontal line at
zero or a main title.
}}
\section{Side Effects}{
a plot is produced on the current graphics device.
}
\seealso{
\code{\link{cox.zph}}, \code{\link{coxph}}
}
\examples{
data(veteran)
vfit <- coxph(Surv(time,status) ~ trt + factor(celltype) + karno +age,
		   data=veteran, x=T)
temp <- cox.zph(vfit)
plot(temp, var=5)      #Look at Karnofsy score
abline(0,0, lty=3)
abline( lm( temp$y[,5] ~ temp$x), lty=4)   #Add the linear fit as well
title(main="VA Lung Study")
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
