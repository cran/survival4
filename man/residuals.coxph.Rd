\name{residuals.coxph}
\title{
Calculate residuals for a coxph fit.
}
\usage{
resid(object,
       type=c("martingale", "deviance", "score", "schoenfeld",
	      "dfbeta", "dfbetas", "scaledsch"),
       collapse)
}
\arguments{
\item{object}{
a coxph object, output from a coxph fit.
}
\item{type}{
character string indicating the type of residual desired;
the default is martingale.
Only enough of the string to determine a unique match is required.
}
\item{collapse}{
Vector indicating which rows to collapse(sum) over.  In time-dependent
models more than one row data can pertain to a single individual.
If there were 4 individuals represented by 3, 1, 2 and 4 rows of data
respectively, then \code{collapse=c(1,1,1, 2, 3,3, 4,4,4,4)} could be used to
obtain per subject rather than per observation residuals.
}}
\value{
For martingale and deviance residuals, the returned object is a vector
with one element for each subject (without \code{collapse}).
For score residuals it is a matrix
with one row per subject and one column per variable.
The row order will match the input data for the original fit.
For Schoenfeld residuals, the returned object is a matrix with one row
for each event and one column per variable.  The rows are ordered by time
within strata, and an attribute \code{strata} is attached that contains the
number of observations in each strata.
The scaled Schoenfeld residuals are used in the cox.zph function.


The score residuals are each individual's contribution to the score vector.
Two transformatons of
this are often more useful: \code{dfbeta} is the approximate change in the
coefficient vector if that observation were dropped,
and 'dfbetas' is the approximate change in the coefficients, scaled by
the standard error for the coefficients.
}
\section{Side Effects}{
For deviance residuals, the status variable may need to be reconstructed.
For score and Shoenfeld residuals, the X matrix will need to be reconstructed.
}
\references{
T.Therneau, P. Grambsch, and T.Fleming. "Martingale based residuals
for survival models", Biometrika, March 1990.
}
\seealso{
\code{\link{coxph}}
}
\examples{
data(jasa1)
fit <- coxph(Surv(start, stop, event) ~ (age + surgery)* transplant,data=jasa1)
mresid <- resid(fit, collapse=jasa1$id)
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
