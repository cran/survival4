\name{pyears}
\title{
Person Years
}
\usage{
pyears(formula, data, weights, subset, na.action, ratetable=survexp.us,
scale=365.25, model=F, x=F, y=F)
}
\arguments{
\item{formula}{
a formula object.  The response variable will be a vector of follow-up times
for each subject, or a Surv object containing the survival time and an
event indicator.
The predictors consist of optional grouping variables
separated by + operators (exactly as in \code{survfit}), time-dependent grouping
variables such as age (specified with \code{tcut}), and optionally a
\code{ratetable()} term.  This latter matches each subject to his/her expected
cohort.
}
\item{data,}{
as in other modeling routines.
Weights are case weights.
}
\item{ratetable}{
a table of event rates, such as survexp.uswhite.
}
\item{scale}{
a scaling for the results.  As most rate tables are in units/day, the
default value of 365.25 causes the output to be reported in years.
}
\item{model,}{
flags to control what is returned.  If any of these is true, then the
model frame, the model matrix, and/or the vector of response times will be
returned as components of the final result, with the same names as the
flag arguments.
}}
\value{
a list with components

\item{pyears}{
an array containing the person-years of exposure. (Or other units, depending
on the rate table and the scale).
}
\item{n}{
an array containing the number of subjects who contribute time to each cell
of the pyears array.
}
\item{event}{
an array containing the observed number of events.  This will be present only
if the resonse variable is a Surv object.
}
\item{expeced}{
an array containing the expected number of events.  This will be present only
if there was a ratetable term.
}
\item{offtable}{
the number of person-years of exposure in the cohort that was not part of
any cell in the pyears array.  This is often useful as an error check; if
there is a mismatch of units between two variables, nearly all the person
years may be off table.
}
\item{call}{
an image of the call to the function.
}
\item{na.action}{
the na.action attribute contributed by an na.action routine, if any.
}}
\description{
Because pyears may have several time variables, it is necessary that all
of them be in the same units.  For instance in the call
 py <- pyears(futime ~ rx + ratetable(age=age, sex=sex, year=entry.dt))
with a ratetable whose natural unit is days, it is important that futime,
age and entry.dt all be in days.  Given the wide range of possible inputs,
it is difficult for the routine to do sanity checks of this aspect.
A special function \code{tcut} is needed to specify time-dependent cutpoints.
For instance, assume that age is in years, and that the desired final
arrays have as one of their margins the age groups 0-2, 2-10, 10-25, and 25+.
A subject who enters the study at age 4 and remains under observation for
10 years will contribute follow-up time to both the 2-10 and 10-25
subsets.  If \code{cut(age, c(0,2,10,25,100))} were used in the formula, the
subject would be classifed according to his starting age only.  The tcut
function has the same arguments as cut, but produces a different output
object which allows the pyears function to correctly track the subject.
The results of pyears() are normally used as input to further calculations.
The print routine, therefore, is designed to give only a summary of the
table.
}
\seealso{
\code{\link{ratetable}}, \code{\link{survexp}}, \code{\link{Surv}}
}
\keyword{survival}
% Converted by Sd2Rd version 0.2-a3.
