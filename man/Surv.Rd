\name{Surv}
\title{
Package a survival variable
}
\usage{
Surv(time, event)  or Surv(time, time2, event)
}
\arguments{
\item{time}{
for right censored data, this is the follow up time.  For interval data, the
first argument is the starting time for the interval.
}
\item{event}{
The status indicator, normally 0=alive, 1=dead.  Other choices are T/F
(TRUE = death) or 1/2 (2=death).
For interval censored data, the status indicator is 0=right censored,
1= event at \code{time}, 2=left censored, 3=interval censored.
}
\item{time2}{
For interval censored  or counting process data only, the ending time of the interval.
Intervals are
assummed to be open on the left and closed on the right, (start, end].
For counting process data,
\code{event} marks whether an event occured at the end of the interval.
}
\item{type}{
one of left, right, counting, interval, or interval2.
If this is not specified, the
default is either right or counting, depending on whether the \code{time2}
argument is absent or present, respectively.
}
\item{origin}{
for counting process data, the hazard function origin.  This is most often
used in conjunction with a model containing time dependent strata in order
to align the subjects properly when they cross over from one strata to
another.
}}
\value{
An object of class 'Surv'.  There are methods for \code{print}, \code{is.na}, and
subscripting survival objects.  To include a survival object inside a
data frame, use the \code{I()} function.  Surv objects are implimented as
a matrix of 2 or 3 columns.
}
\section{METHOD}{
In theory it is possible to represent interval censored data without a
third column containing the explicit status.  Exact, right censored,
left censored and interval censored observation would be represented as
intervals of (a,a), (a, infinity), (-infinity,b), and (a,b) respectively;
each specifing the interval within which the event is known to have occured.
Infinity is, of course, impractical in a computer routine.
If \code{type} is "interval2" then the representation given above is assumed,
with NA taking the place of infinity.  If type='interval' then an explicit
status code must be given in the third argument.
If the status code is 0, 1 or 2, then the relevant information is assumed to
be contained in \code{time},  the value in \code{time2} is ignored, and the second
column of the result will contain a placeholder.
At present, all of the methods that handle interval censored data are
parametric models, so the distinction between open and closed intervals
is unimportant.  The distinction is important for counting process data and
the Cox model.
}
\examples{
data(aml)
Surv(aml$time, aml$status)
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
