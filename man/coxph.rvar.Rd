\name{coxph.rvar}
\title{
Robust variance for a Cox model
}
\usage{
coxph.rvar(fit, collapse)
}
\arguments{
\item{fit}{
a coxph object, i.e., the result of fitting a Cox model.
}
\item{collapse}{
if the original data contained correlated observations, e.g., multiple
data rows per subject, then this argument contains the id vector that
identifies the subgroups.
}
\value{
a copy of the input, with two components added
}
\item{robust.var}{
the robust variance estimate.
}
\item{rcall}{
the call to this function.
}}
\section{Side Effects}{
the print and summary methods for coxph recognize and use the robust
variance. The global likelihood ratio and score statistics are
uneffected, but the global Wald test will now be based on the robust
estimator.
}
\description{
Let r be the matrix of infinitesimal influence functions, i.e.,
r <- residuals(fit, type='dbeta').  Then the robust variance is
v <- t(r) \%*\% r.  If there are correlated observations, the appropriate rows
or r are first summed, and v is based on the reduced r matrix.  There is
an obvious connection with the ordinary and group jackknife estimates.
}
\seealso{
\code{\link{coxph}}
}
\examples{
data(ovarian)
fit <- coxph(Surv(futime, fustat) ~ age + rx +ecog.ps, data=ovarian)
fit2 <- coxph.rvar(fit)
summary(fit2)
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
