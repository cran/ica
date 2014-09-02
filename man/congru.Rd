\name{congru}
\alias{congru}
\title{
Tucker's congruence coefficient
}
\description{
Tucker's congruence coefficient is the unadjusted (or uncentered) correlation coefficient between two vectors \code{x} and \code{y}. Ranges between -1 (perfect negative association) and 1 (perfect positive association) with a value of 0 denoting no linear relationship.
}
\usage{
congru(x,y)
}
\arguments{
  \item{x}{
Vector of dimension \eqn{n \times 1}.
}
  \item{y}{
Vector of dimension \eqn{n \times 1}.
}
}
\value{
Returns a scalar (the congruence coefficient). 
}
\references{
Tucker, L.R. (1951). \emph{A method for synthesis of factor analysis studies} (Personnel Research Section Report No. 984). Washington, DC: Department of the Army.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
Tucker's congruence coefficient is defined as \deqn{r = \frac{\sum_{i=1}^{n}x_{i}y_{i}}{ \sqrt{\sum_{i=1}^{n}x_{i}^{2}\sum_{i=1}^{n}y_{i}^{2}} } }
where \eqn{x_{i}} and \eqn{y_{i}} denote the \eqn{i}-th elements of \code{x} and \code{y}.
}
\note{
If inputs \code{x} and \code{y} are matrices, inputs will be vectorized before calculating \emph{r}.
}
\examples{
##########   EXAMPLE   ##########

A<-rnorm(100)
B<-rnorm(100)
C<-A*5
D<-A*(-0.5)
congru(A,B)
congru(A,C)
congru(A,D)

}
