\name{is.ratetable}
\title{Is an Object of Class ``ratetable''?}
\usage{
is.ratetable(x, verbose=F)
}
\arguments{
  \item{x}{the object to be verified.}
  \item{verbose}{
    If \code{TRUE} and the object is not a ratetable, then return a
    character string describing the way(s) in which \code{x} fails to be
    a proper ratetable object.}
}
\description{
  Verify that an object is of class \code{ratetable}.  The function
  verifies not only the \code{class} attribute, but the structure of the
  object.
}
\value{
  Returns \code{TRUE} if \code{x} is a ratetable, and \code{FALSE} or a
  description if it is not.  Rate tables are used by the \code{pyears}
  and \code{survexp} functions, and normally contain death rates for
  some population, categorized by age, sex, or other variables.  They
  have a fairly rigid structure, and the \code{verbose} option can help
  in creating a new rate table. 
}
\seealso{
  \code{\link{pyears}}, \code{\link{survexp}}
}
