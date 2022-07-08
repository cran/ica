\name{icaimax}
\alias{icaimax}
\title{
ICA via Infomax Algorithm
}
\description{
Computes ICA decomposition using Bell and Sejnowski's (1995) Information-Maximization (Infomax) approach with various options.
}
\usage{
icaimax(X, nc, center = TRUE, maxit = 100, tol = 1e-6, Rmat = diag(nc), 
        alg = "newton", fun = "tanh", signs = rep(1, nc), signswitch = TRUE, 
        rate = 1, rateanneal = NULL)
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
Algorithm to use: \code{alg="newton"} for Newton iteration, and \code{alg="gradient"} for gradient descent.
}
  \item{fun}{
Nonlinear (squashing) function to use for algorithm: \code{fun="tanh"} for hyperbolic tangent, \code{fun="log"} for logistic, and \code{fun="ext"} for extended Infomax.  
}
  \item{signs}{
Vector of length \code{nc} such that \code{signs[j] = 1} if j-th component is super-Gaussian and \code{signs[j] = -1} if j-th component is sub-Gaussian. Only used if \code{fun="ext"}. Ignored if \code{signswitch=TRUE}.
}
  \item{signswitch}{
If \code{TRUE}, the \code{signs} vector is automatically determined from the data; otherwise a confirmatory ICA decomposition is calculated using input \code{signs} vector. Only used if \code{fun="ext"}.
}
  \item{rate}{
Learing rate for gradient descent algorithm. Ignored if \code{alg="newton"}.  
}
  \item{rateanneal}{
Annealing angle and proportion for gradient descent learing rate (see Details). Ignored if \code{alg="newton"}.  
}


}
\value{
\item{S}{Matrix of source signal estimates (\code{S = Y \%*\% R}).}
\item{M}{Estimated mixing matrix.}
\item{W}{Estimated unmixing matrix (\code{W = crossprod(R, Q)}).}
\item{Y}{Whitened data matrix.}
\item{Q}{Whitening matrix.}
\item{R}{Orthogonal rotation matrix.}
\item{vafs}{Variance-accounted-for by each component.}
\item{iter}{Number of algorithm iterations.}
\item{alg}{Algorithm used (same as input).}
\item{fun}{Contrast function (same as input).}
\item{signs}{Component signs (same as input).}
\item{rate}{Learning rate (same as input).}
\item{converged}{Logical indicating if algorithm converged.}
}
\references{
Bell, A.J. & Sejnowski, T.J. (1995). An information-maximization approach to blind separation and blind deconvolution. \emph{Neural Computation, 7}(6), 1129-1159. \doi{10.1162/neco.1995.7.6.1129}

Helwig, N.E. & Hong, S. (2013). A critique of Tensor Probabilistic Independent Component Analysis: Implications and recommendations for multi-subject fMRI data analysis. \emph{Journal of Neuroscience Methods, 213}(2), 263-273. \doi{https://doi.org/10.1016/j.jneumeth.2012.12.009}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
\bold{ICA Model}
The ICA model can be written as \code{X = tcrossprod(S, M) + E}, where \code{S} contains the source signals, \code{M} is the mixing matrix, and \code{E} contains the noise signals. Columns of \code{X} are assumed to have zero mean. The goal is to find the unmixing matrix \code{W} such that columns of \code{S = tcrossprod(X, W)} are independent as possible.

\bold{Whitening}
Without loss of generality, we can write \code{M = P \%*\% R} where \code{P} is a tall matrix and \code{R} is an orthogonal rotation matrix. Letting \code{Q} denote the pseudoinverse of \code{P}, we can whiten the data using \code{Y = tcrossprod(X, Q)}. The goal is to find the orthongal rotation matrix \code{R} such that the source signal estimates \code{S = Y \%*\% R} are as independent as possible. Note that \code{W = crossprod(R, Q)}.

\bold{Infomax}
The Infomax approach finds the orthogonal rotation matrix \code{R} that (approximately) maximizes the joint entropy of a nonlinear function of the estimated source signals. See Bell and Sejnowski (1995) and Helwig and Hong (2013) for specifics of algorithms.

}

\seealso{
\code{\link{icafast}} for FastICA

\code{\link{icajade}} for ICA via JADE
}

\examples{
##########   EXAMPLE 1   ##########

# generate noiseless data (p == r)
set.seed(123)
nobs <- 1000
Amat <- cbind(icasamp("a", "rnd", nobs), icasamp("b", "rnd", nobs))
Bmat <- matrix(2 * runif(4), nrow = 2, ncol = 2)
Xmat <- tcrossprod(Amat, Bmat)

# ICA via Infomax with 2 components
imod <- icaimax(Xmat, nc = 2)
acy(Bmat, imod$M)
cor(Amat, imod$S)



##########   EXAMPLE 2   ##########

# generate noiseless data (p != r)
set.seed(123)
nobs <- 1000
Amat <- cbind(icasamp("a", "rnd", nobs), icasamp("b", "rnd", nobs))
Bmat <- matrix(2 * runif(200), nrow = 100, ncol = 2)
Xmat <- tcrossprod(Amat, Bmat)

# ICA via Infomax with 2 components
imod <- icaimax(Xmat, nc = 2)
cor(Amat, imod$S)



##########   EXAMPLE 3   ##########

# generate noisy data (p != r)
set.seed(123)
nobs <- 1000
Amat <- cbind(icasamp("a", "rnd", nobs), icasamp("b", "rnd", nobs))
Bmat <- matrix(2 * runif(200), 100, 2)
Emat <- matrix(rnorm(10^5), nrow = 1000, ncol = 100)
Xmat <- tcrossprod(Amat,Bmat) + Emat

# ICA via Infomax with 2 components
imod <- icaimax(Xmat, nc = 2)
cor(Amat, imod$S)

}
