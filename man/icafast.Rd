\name{icafast}
\alias{icafast}
\title{
ICA via FastICA Algorithm
}
\description{
Computes ICA decomposition using Hyvarinen's (1999) FastICA algorithm with various options.
}
\usage{
icafast(X, nc, center = TRUE, maxit = 100, tol = 1e-6,
        Rmat = diag(nc), alg = c("par", "def"),
        fun = c("logcosh", "exp", "kur"), alpha = 1)
}
\arguments{
  \item{X}{
Data matrix with \code{n} rows (samples) and \code{p} columns (variables).
}
  \item{nc}{
Number of components to extract.
}
  \item{center}{
If \code{TRUE}, columns of \code{X} are mean-centered before ICA decomposition.
}
  \item{maxit}{
Maximum number of algorithm iterations to allow.
}
  \item{tol}{
Convergence tolerance.
}
  \item{Rmat}{
Initial estimate of the \code{nc}-by-\code{nc} orthogonal rotation matrix.
}
  \item{alg}{
Algorithm to use: \code{alg="par"} to estimate all \code{nc} components in parallel (default) or \code{alg="def"} for deflation estimation (i.e., projection pursuit).
}
  \item{fun}{
Contrast function to use for negentropy approximation.  
}
  \item{alpha}{
Tuning parameter for "logcosh" contrast function (1 <= \code{alpha} <= 2).
}
}
\value{
\item{S}{Matrix of source signal estimates (\code{S=Y\%*\%R}).}
\item{M}{Estimated mixing matrix.}
\item{W}{Estimated unmixing matrix (\code{W=crossprod(R,Q)}).}
\item{Y}{Whitened data matrix.}
\item{Q}{Whitening matrix.}
\item{R}{Orthogonal rotation matrix.}
\item{vafs}{Variance-accounted-for by each component.}
\item{iter}{Number of algorithm iterations.}
\item{alg}{Algorithm used (same as input).}
\item{fun}{Contrast function (same as input).}
\item{alpha}{Tuning parameter (same as input).}
}
\references{
Helwig, N.E. & Hong, S. (2013). A critique of Tensor Probabilistic Independent Component Analysis: Implications and recommendations for multi-subject fMRI data analysis. \emph{Journal of Neuroscience Methods, 213}, 263-273.

Hyvarinen, A. (1999). Fast and robust fixed-point algorithms for independent component analysis. \emph{IEEE Transactions on Neural Networks, 10}, 626-634.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
\bold{ICA Model}
The ICA model can be written as \code{X=tcrossprod(S,M)+E}, where columns of \code{S} contain the source signals, \code{M} is the mixing matrix, and columns of \code{E} contain the noise signals. Columns of \code{X} are assumed to have zero mean. The goal is to find the unmixing matrix \code{W} such that columns of \code{S=tcrossprod(X,W)} are independent as possible.

\bold{Whitening}
Without loss of generality, we can write \code{M=P\%*\%R} where \code{P} is a tall matrix and \code{R} is an orthogonal rotation matrix. Letting \code{Q} denote the pseudoinverse of \code{P}, we can whiten the data using \code{Y=tcrossprod(X,Q)}. The goal is to find the orthongal rotation matrix \code{R} such that the source signal estimates \code{S=Y\%*\%R} are as independent as possible. Note that \code{W=crossprod(R,Q)}.

\bold{FastICA}
The FastICA algorithm finds the orthogonal rotation matrix \code{R} that (approximately) maximizes the negentropy of the estimated source signals. Negentropy is approximated using \deqn{J(s) = [E\{G(s)\}-E\{G(z)\} ]^2} where \emph{E} denotes the expectation, \emph{G} is the contrast function, and \emph{z} is a standard normal variable. See Hyvarinen (1999) for specifics of fixed-point algorithm.

}

\examples{
##########   EXAMPLE 1   ##########

# generate noiseless data (p==r)
set.seed(123)
nobs <- 1000
Amat <- cbind(icasamp("a","rnd",nobs),icasamp("b","rnd",nobs))
Bmat <- matrix(2*runif(4),2,2)
Xmat <- tcrossprod(Amat,Bmat)

# ICA via FastICA with 2 components
imod <- icafast(Xmat,2)
acy(Bmat,imod$M)
cor(Amat,imod$S)



##########   EXAMPLE 2   ##########

# generate noiseless data (p!=r)
set.seed(123)
nobs <- 1000
Amat <- cbind(icasamp("a","rnd",nobs),icasamp("b","rnd",nobs))
Bmat <- matrix(2*runif(200),100,2)
Xmat <- tcrossprod(Amat,Bmat)

# ICA via FastICA with 2 components
imod <- icafast(Xmat,2)
cor(Amat,imod$S)



##########   EXAMPLE 3   ##########

# generate noisy data (p!=r)
set.seed(123)
nobs <- 1000
Amat <- cbind(icasamp("a","rnd",nobs),icasamp("b","rnd",nobs))
Bmat <- matrix(2*runif(200),100,2)
Emat <- matrix(rnorm(10^5),1000,100)
Xmat <- tcrossprod(Amat,Bmat)+Emat

# ICA via FastICA with 2 components
imod <- icafast(Xmat,2)
cor(Amat,imod$S)

}
