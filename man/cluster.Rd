\name{cluster}
\title{Identify clusters}
\usage{
cluster(x)
}
\arguments{
\item{x}{A character, factor, or numeric variable.}
}
\description{
  This is a special function used in the context of the Cox model.  It
  identifies correlated groups of observations, and is used on the right
  hand side of a formula.
}
\value{
  The function's only action is semantic, to mark a variable as the
  cluster indicator.
}
\seealso{
  \code{\link{coxph}}
}
\examples{
data(ovarian)
coxph(Surv(futime, fustat) ~ age + cluster(group))
}
