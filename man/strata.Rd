\name{strata}
\title{
Identify strata variables.
}
\usage{
strata(\dots, na.group=F)
}
\arguments{
\item{...}{
Any number of variables.  All must be the same length.
}
\item{na.group}{
if set to T, then missing values are treated as a distinct level of each
variable.
}}
\value{
a new factor, whose levels are all possible combinations of the factors
supplied as arguments.
}
\description{
The result is identical to the interaction() function, but for the
labeling of the factors (strata is more verbose).
}
\seealso{
\code{\link{coxph}}
}
\examples{
data(ovarian)
coxph(Surv(futime, fustat) ~ age + strata(rx))
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
