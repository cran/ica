\name{acyidx}
\alias{acyidx}
\title{
Amari-Cichocki-Yang (ACY) performance index
}
\description{
The ACY index is a measure of dissimilarity between two (nonsingular) matrices \code{A} and \code{B} of the same dimension. The ACY index: (a) is invariant to permutation and rescaling of the columns of \code{A} and \code{B}, (b) ranges between 0 and n-1, and (c) equals 0 if and only if \code{A} and \code{B} are identical up to column permutations and rescalings.
}
\usage{
acyidx(X,Y)
}
\arguments{
  \item{X}{
Nonsingular matrix of dimension \eqn{n \times n}.
}
  \item{Y}{
Nonsingular matrix of dimension \eqn{n \times n}.
}
}
\value{
Returns a scalar (the ACY index). 
}
\references{
Amari, S., Cichocki, A., & Yang, H.H. (1996). A new learning algorithm for blind signal separation. In D. S. Touretzky, M. C. Mozer, and M. E. Hasselmo (Eds.), \emph{Advances in Neural Information Processing Systems, 8}. Cambridge, MA: MIT Press.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
The acyidx index is defined as \deqn{\frac{1}{2n}\sum_{i=1}^{n}\left(\frac{\sum_{j=1}^{n}|a_{ij}|}{\max_{j}|a_{ij}|}-1\right) + \frac{1}{2n}\sum_{j=1}^{n}\left(\frac{\sum_{i=1}^{n}|a_{ij}|}{\max_{i}|a_{ij}|}-1\right)  }
where \eqn{a_{ij} = (\mathbf{Y}^{-1}\mathbf{X})_{ij}}. 
}
\section{Warnings }{
If \code{B} is singular, function will produce an error.
}

\examples{
##########   EXAMPLE   ##########

set.seed(123)
A<-matrix(runif(16),4,4)
B<-matrix(runif(16),4,4)
C<-A[,c(3,1,2,4)]\%*\%diag(1:4)
acyidx(A,B)
acyidx(A,C)

}
