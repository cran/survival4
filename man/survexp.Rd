\name{survexp}
\title{
Compute Expected Survival
}
\usage{
survexp(formula, data, weights, subset, na.action,
 times, cohort=T, conditional=T,
 ratetable=survexp.us, scale=1, se.fit, model=F, x=F, y=F)
}
\arguments{
\item{formula}{
a formula object.  The response variable will be a vector of follow-up times,
and is optional.  The predictors will consist of optional grouping variables
separated by + operators (exactly as in \code{survfit}), along with a
\code{ratetable()} term.  This latter matches each subject to his/her expected
cohort.
}
\item{data,}{
as in other modeling routines.
Weights are currently ignored.
}
\item{times}{
an optional vector of times at which the resulting survival curve should
be evaluated.  If absent, the result will be reported for each unique value
of the vector of follow-up times.
}
\item{cohort}{
If false, each subject is treated as a subgroup of size 1.
}
\item{conditional}{
If \code{y} is missing in the formula, this argument is ignored.  Otherwise it
is an indicator of whether y includes death times, which leads to conditional
expected survival, or y includes only the potential censoring times.
}
\item{ratetable}{
a table of event rates, such as survexp.uswhite, or a fitted Cox model.
}
\item{scale}{
a scaling for the results.  As most rate tables are in units/day, a
value of 365.24 would cause the output to be reported in years.
}
\item{npoints}{
calculate intermediate results at npoints values, evenly spaced on the range
of \code{y}.  The usual (exact) calculation is done at each unique 'y' value;
for very large data sets this may incur too much storage for the scratch
array.
For a prediction from a Cox model this arument is ignored.
}
\item{se.fit}{
compute the standard error of the predicted survival.
The default is to compute this whenever the routine can, which at this time
is only for the Ederer method and a Cox model as the rate table.
}
\item{model,}{
flags to control what is returned.  If any of these is true, then the
model frame, the model matrix, and/or the vector of response times will be
returned as components of the final result, with the same names as the
flag arguments.
}}
\value{
if cohort=T an object of class \code{survexp}, otherwise a vector of per-subject
expected survival values.  The former contains the number of subjects at
risk and the expected survival for the cohort at each requested time.
}
\description{
Individual expected survival is ususally used in models or testing, to
\code{correct} for the age and sex composition of a group of subjects.  For
instance, assume that birth date, entry date onto the study,sex and
actual survival time are all known for a group of subjects.
The uswhite population tables contain expected death rates
based on calendar year, sex and age.  Then
haz <- -log(survexp(death.time ~ ratetable(sex=sex, year=entry.dt, age=(birth.dt-entry.dt)), cohort=F))
gives for each subject the total hazard experienced up to their observed
death time or censoring time.
This probability can be used as a rescaled time value in models:
glm(status ~ 1 + offset(log(haz)), family=poisson)
glm(status ~ x + offset(log(haz)), family=poisson)
In the first model, a test for intercept=0 is the one sample log-rank
test of whether the observed group of subjects has equivalent survival to
the baseline population.  The second model tests for an effect of variable
\code{x} after adjustment for age and sex.


Cohort survival is used to produce an overall survival curve.  This is then
added to the Kaplan-Meier plot of the study group for visual comparison
between these subjects and the population at large.  There are three common
methods of computing cohort survival.
In the "exact method" of Ederer the cohort is not censored; this corresponds
to having no response variable in the formula.  Hakulinen recommends censoring
the cohort at the anticipated censoring time of each patient, and Verhuel
recommends censoring the cohort at the actual observation time of each
patient.
The last of these is the conditional method.
These are obtained by using the respective time values as the
follow-up time or response in the formula.
}
\references{
G. Berry.  The analysis of mortality by the subject-years method.
Biometrics 1983, 39:173-84.
F Ederer, L Axtell, and S Cutler.  The relative survival rate: a statistical
methodology. Natl Cnacer Inst Monogr 1961, 6:101-21.
T. Hakulinen.  Cancer survival corrected for heterogeneity in patient
withdrawal.  Biometrics 1892, 38:933.
H. Verheul, E. Dekker, P. Bossuyt, A. Moulijn, and A. Dunning.  Backround
mortality in clinical survival studies.  Lancet 1993, 341:872-5.
}
\seealso{
\code{\link{survfit}}, \code{\link{survexp.us}}, \code{\link{survexp.fit}}, \code{\link{personyr}}, \code{\link{date}}
}
\examples{
\dontrun{
library(ratetables)
data(jasa)
efit <- survexp( ~ ratetable(sex=sex, year=entry.dt, age=entry.dt-birth.dt))
plot(survfit(Surv(futime, status) ~1))
lines(efit)
}
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
