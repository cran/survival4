\name{cox.zph}
\title{
Test the proportional hazards assumption
}
\usage{
cox.zph(fit, transform='km', global=T)
}
\arguments{
\item{fit}{
the result of fitting a Cox regression model, using the coxph function.
}
\item{transform}{
survival times transformed before the test is performed.  Possible values
are the character strings 'km', 'rank', 'identity', or a function of one argument.
}
\item{global}{
should a global chisquare test be done, in additoin to the
 per-variable tests.
}
\item{x}{
if true, then the result will be a list containing the test table (a matrix),
x and y.  If false then only the test table is returned.
}
\value{
an object of class \code{cox.zph}, with components:
}
\item{table}{
a matrix with one row for each variable, and optionally a last row for
the global test.
Columns of the matrix contain the correlation coefficient between transformed
survival time and the scaled Schoenfeld residuals, a chisquare,
and the two-sided p value.
For the global test there is no appropriate correlation, so an NA is
entered into the matrix as a placeholder.
}
\item{x}{
the transformed time axis.
}
\item{y}{
the matrix of scaled Schoenfeld residuals.  There will be one column per
variable and one row per event.  The row labels contain the original event
times, which will be the same as \code{x} for the identity transform.
}
\item{call}{
the calling sequence for the routine.
The computations require the original \code{x} matrix of the Cox model fit.
Thus it saves time if the \code{x=T} option is used in coxph.
This function would usually be executed in conjunction with a plot of the
scaled Schoenfeld residuals versus survival time.  The plot gives an
estimate of the time-dependent coefficient beta(t).
}}
\references{
Terry Therneau, author of local function.
}
\note{The \code{plot} method requires the \code{splines} library}
\examples{
data(ovarian)
fit <- coxph( Surv(futime, fustat) ~ age + surgery, jasa, x=T)
temp<- cox.zph(fit, x=T)
print(temp)                  #display the results
plot(temp)                   #plot curves
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
