\name{Rank distance correlation}
\alias{rdcor}
\title{
Rank distance correlation
}

\description{
Rank distance correlation.
}

\usage{
rdcor(y, x)
}

\arguments{
\item{y}{
A numerical vector.
}
\item{x}{
A numerical vector or a numerical matrix.
}
}

\details{
This computes the rank distance correlation between y and x, or between y and each column of x.
}

\value{
A vector with the rank distance correlation().
}

\references{
Shi H., Drton M. and Han F. (2022). Distribution-free consistent independence tests via center-outward ranks and signs. Journal of the American Statistical Association, 117(537): 395--410.

Zhang Q. (2025). On the connections between Chatterjee's correlation and rank distance correlation.
Journal of Nonparametric Statistics, 1--18.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{ \link{rdcor.test}, \link{rdcor.mat}
}
}

\examples{
y <- iris[, 1]
x <- matrix( rnorm(150 * 10), ncol = 10 )
rdcor(y, x)
}

