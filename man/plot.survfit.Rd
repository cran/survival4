\name{plot.survfit}
\title{
Plot method for survfit.
}
\usage{
plot.survfit(survfit, conf.int=<<see below>>, mark.time=T,
 mark=3, col=1, lty=1, lwd=1, cex=1, log=F, yscale=1, xscale=1, xlab="",
 ylab="", xaxs='i', \dots)
}
\arguments{
\item{survfit}{
structure returned by survfit.
}
\item{conf.int}{
determines whether confidence intervals will be plotted.  The default is to
do so if there is only 1 curve, i.e., no strata.
}
\item{mark.time}{
controls the labeling of the curves.  If set to False, no labeling is done.
If True, then curves are marked at each censoring time.  If mark.time is a
numeric vector, then curves are marked at the specified time points.
}
\item{mark}{
vector of mark parameters, which will be used to label the curves.
The \code{lines} help file contains examples of the possible marks.
The vector is reused cyclically if it is shorter than the number of curves.
}
\item{col}{
vector of colors.  The default value is 1.
}
\item{lty}{
vector of line types. The default value is 1.
}
\item{lwd}{
vector of line widths. The default value is 1.
}
\item{cex}{
parameter available to change the size of "mark".
Not a vector; all marks have the same size.
}
\item{log}{
logical value: should the y axis be on a log scale?
}
\item{yscale}{
will be used to multiply the labels on the y axis.
A value of 100, for instance, would be used to give a percent scale.
Only the labels are
changed, not the actual plot coordinates, so that adding a curve with
"lines(surv.exp(\dots{}))", say, will perform as it did without the yscale arg.
}
\item{yscale}{
will be used in a similar manner for labels on the x axis.  A value of
365.25 will give labels in years instead of the original days.
}
\item{xlab}{
label given to the x-axis.
}
\item{ylab}{
label given to the y-axis.
}
\item{xaxs}{
the x axis style, as listed in \code{par}.  The S default option of "r" leads to
curves that are "too far" from the y axis.  This is, of course, just a matter
of esthetic opinion.
}}
\value{
a list with components x and y, containing the coordinates of the last point
on each of the curves.  This may be useful for labeling.
}
\section{Side Effects}{
A plot of survival curves is produced, one curve for each strata.
}
\examples{
data(ovarian)
plot.survfit(survfit(Surv(futime,fustat)~rx,data=ovarian))
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.


