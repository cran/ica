\name{acy}
\alias{acy}
\title{
Amari-Cichocki-Yang Error
}
\description{
The Amari-Cichocki-Yang (ACY) error is an asymmetric measure of dissimilarity between two nonsingular matrices \code{X} and \code{Y}. The ACY error: (a) is invariant to permutation and rescaling of the columns of \code{X} and \code{Y}, (b) ranges between 0 and \code{n-1}, and (c) equals 0 if and only if \code{X} and \code{Y} are identical up to column permutations and rescalings.
}
\usage{
acy(X,Y)
}
\arguments{
  \item{X}{
Nonsingular matrix of dimension \eqn{n \times n} (test matrix).
}
  \item{Y}{
Nonsingular matrix of dimension \eqn{n \times n} (target matrix).
}
}
\value{
Returns a scalar (the ACY error).
}
\references{
Amari, S., Cichocki, A., & Yang, H.H. (1996). A new learning algorithm for blind signal separation. In D. S. Touretzky, M. C. Mozer, and M. E. Hasselmo (Eds.), \emph{Advances in Neural Information Processing Systems, 8}. Cambridge, MA: MIT Press.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
The ACY error is defined as \deqn{\frac{1}{2n}\sum_{i=1}^{n}\left(\frac{\sum_{j=1}^{n}|a_{ij}|}{\max_{j}|a_{ij}|}-1\right) + \frac{1}{2n}\sum_{j=1}^{n}\left(\frac{\sum_{i=1}^{n}|a_{ij}|}{\max_{i}|a_{ij}|}-1\right)  }
where \eqn{a_{ij} = (\mathbf{Y}^{-1}\mathbf{X})_{ij}}. 
}
\section{Warnings }{
If \code{Y} is singular, function will produce an error.
}

\examples{
##########   EXAMPLE   ##########

set.seed(1)
X <- matrix(runif(16),4,4)
Y <- matrix(runif(16),4,4)
Z <- X[,c(3,1,2,4)]\%*\%diag(1:4)
acy(X,Y)
acy(X,Z)

}
