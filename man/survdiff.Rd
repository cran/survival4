\name{survdiff}
\title{
Test Survival Curve Differences
}
\usage{
survdiff(formula, data,  rho=0, subset)
}
\arguments{
\item{formula}{
a formula expression as for other survival models, of the form
\code{Surv(time, status) ~ predictors}.  For a one-sample test, the predictors
must consist of a single \code{offset(sp)} term, where sp is a vector giving the
survival probability of each subject.  For a k-sample test, each unique
combination of predictors defines a subgroup.
To cause missing values in the predictors to be treated as a separate
group, rather than being omitted, use the \code{strata} function with its
\code{na.group=T} argument.
}
\item{data}{
an optional data frame in which to interpret the variables occurring in the
formula.
}
\item{rho}{
a parameter that controls the type of test.
}
\item{subset}{
subset of the observations to be used in the fit.
}
\value{
a list with components:
}
\item{n}{
the number of subjects in each group.
}
\item{obs}{
the weighted observed number of events in each group.
}
\item{exp}{
the weighted expected number of events in each group.
}
\item{chisq}{
the chisquare statistic for a test of equality.
}}
\section{METHOD}{
This function implements the G-rho family of
Harrington and Fleming (1982), with weights on each death of (S(t))^rho,
where S is the Kaplan-Meier estimate of survival.
When \code{rho = 0} this is the log-rank or Mantel-Haenszel test,
and when \code{rho = 1} it is equivalent to the Peto \& Peto modification
of the Gehan-Wilcoxon test.


If the right hand side of the formula consists only of an offset term,
then a one sample test is done.
}
\references{
Harrington, D. P. and Fleming, T. R. (1982).
A class of rank test procedures for censored survival data.
Biometrika
\bold{69}, 553-566.
}
\seealso{
\code{\link{survdiff.print}}.
}
\examples{
data(ovarian)
survdiff(Surv(futime, fustat) ~ rx)
dontrun{
library(ratetables)
expect <- survexp(entry, birth, sex, futime)
survdiff(Surv(futime, fustat) ~ offset(expect$surv))  #One sample log-rank
}
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
