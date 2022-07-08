\name{icajade}
\alias{icajade}
\title{
ICA via JADE Algorithm
}
\description{
Computes ICA decomposition using Cardoso and Souloumiac's (1993, 1996) Joint Approximate Diagonalization of Eigenmatrices (JADE) approach.
}
\usage{
icajade(X, nc, center = TRUE, maxit = 100, tol = 1e-6, Rmat = diag(nc))
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
\item{converged}{Logical indicating if algorithm converged.}
}
\references{
Cardoso, J.F., & Souloumiac, A. (1993). Blind beamforming for non-Gaussian signals. \emph{IEE Proceedings-F, 140}(6), 362-370. \doi{https://doi.org/10.1049/ip-f-2.1993.0054}
  
Cardoso, J.F., & Souloumiac, A. (1996). Jacobi angles for simultaneous diagonalization. \emph{SIAM Journal on Matrix Analysis and Applications, 17}(1), 161-164. \doi{10.1137/S0895479893259546}

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

\bold{JADE}
The JADE approach finds the orthogonal rotation matrix \code{R} that (approximately) diagonalizes the cumulant array of the source signals. See Cardoso and Souloumiac (1993,1996) and Helwig and Hong (2013) for specifics of the JADE algorithm.
}

\seealso{
\code{\link{icafast}} for FastICA

\code{\link{icaimax}} for ICA via Infomax
}

\examples{
##########   EXAMPLE 1   ##########

# generate noiseless data (p == r)
set.seed(123)
nobs <- 1000
Amat <- cbind(icasamp("a", "rnd", nobs), icasamp("b", "rnd", nobs))
Bmat <- matrix(2 * runif(4), nrow = 2, ncol = 2)
Xmat <- tcrossprod(Amat, Bmat)

# ICA via JADE with 2 components
imod <- icajade(Xmat, nc = 2)
acy(Bmat, imod$M)
cor(Amat, imod$S)



##########   EXAMPLE 2   ##########

# generate noiseless data (p != r)
set.seed(123)
nobs <- 1000
Amat <- cbind(icasamp("a", "rnd", nobs), icasamp("b", "rnd", nobs))
Bmat <- matrix(2 * runif(200), nrow = 100, ncol = 2)
Xmat <- tcrossprod(Amat, Bmat)

# ICA via JADE with 2 components
imod <- icajade(Xmat, nc = 2)
cor(Amat, imod$S)



##########   EXAMPLE 3   ##########

# generate noisy data (p != r)
set.seed(123)
nobs <- 1000
Amat <- cbind(icasamp("a", "rnd", nobs), icasamp("b", "rnd", nobs))
Bmat <- matrix(2 * runif(200), 100, 2)
Emat <- matrix(rnorm(10^5), nrow = 1000, ncol = 100)
Xmat <- tcrossprod(Amat,Bmat) + Emat

# ICA via JADE with 2 components
imod <- icajade(Xmat, nc = 2)
cor(Amat, imod$S)

}
