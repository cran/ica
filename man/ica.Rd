\name{ica}
\alias{ica}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ICA via FastICA, Infomax, or JADE
}
\description{
Computes ICA decomposition using Hyvarinen's (1999) FastICA algorithm, Bell and Sejnowski's (1995) Information-Maximization (Infomax) algorithm, or Cardoso and Souloumiac's (1993, 1996) Joint Approximate Diagonalization of Eigenmatrices (JADE) algorithm.
}
\usage{
ica(X, nc, method = c("fast", "imax", "jade"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data matrix with \code{n} rows (samples) and \code{p} columns (variables).
}
  \item{nc}{
Number of components to extract.  
}
  \item{method}{
Method for decomposition.  
}
  \item{...}{
Additional arguments to be passed to other ICA functions (see Details).  
}
}
\details{
\bold{ICA Model}
The ICA model can be written as \code{X = tcrossprod(S, M) + E}, where \code{S} contains the source signals, \code{M} is the mixing matrix, and \code{E} contains the noise signals. Columns of \code{X} are assumed to have zero mean. The goal is to find the unmixing matrix \code{W} such that columns of \code{S = tcrossprod(X, W)} are independent as possible.

\bold{Whitening}
Without loss of generality, we can write \code{M = P \%*\% R} where \code{P} is a tall matrix and \code{R} is an orthogonal rotation matrix. Letting \code{Q} denote the pseudoinverse of \code{P}, we can whiten the data using \code{Y = tcrossprod(X, Q)}. The goal is to find the orthongal rotation matrix \code{R} such that the source signal estimates \code{S = Y \%*\% R} are as independent as possible. Note that \code{W = crossprod(R, Q)}.

\bold{Method}
This is a wrapper function for the functions \code{\link{icafast}}, \code{\link{icaimax}}, or \code{\link{icajade}}. See the corresponding function for details on the method, as well as the available arguments (handled by the \code{...} argument).
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
\item{...}{Other arguments (if \code{method = "fast"} or \code{method = "imax"}).}
}
\references{
Bell, A.J. & Sejnowski, T.J. (1995). An information-maximization approach to blind separation and blind deconvolution. \emph{Neural Computation, 7}(6), 1129-1159. \doi{10.1162/neco.1995.7.6.1129}

Cardoso, J.F., & Souloumiac, A. (1993). Blind beamforming for non-Gaussian signals. \emph{IEE Proceedings-F, 140}(6), 362-370. \doi{https://doi.org/10.1049/ip-f-2.1993.0054}
  
Cardoso, J.F., & Souloumiac, A. (1996). Jacobi angles for simultaneous diagonalization. \emph{SIAM Journal on Matrix Analysis and Applications, 17}(1), 161-164. \doi{10.1137/S0895479893259546}

Helwig, N.E. & Hong, S. (2013). A critique of Tensor Probabilistic Independent Component Analysis: Implications and recommendations for multi-subject fMRI data analysis. \emph{Journal of Neuroscience Methods, 213}(2), 263-273. \doi{https://doi.org/10.1016/j.jneumeth.2012.12.009}

Hyvarinen, A. (1999). Fast and robust fixed-point algorithms for independent component analysis. \emph{IEEE Transactions on Neural Networks, 10}(3), 626-634. \doi{10.1109/72.761722}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{icafast}} for ICA via FastICA

\code{\link{icaimax}} for ICA via Infomax

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

# ICA via different algorithms
imod.fast <- ica(Xmat, nc = 2)
imod.imax <- ica(Xmat, nc = 2, method = "imax")
imod.jade <- ica(Xmat, nc = 2, method = "jade")

# compare mixing matrix recovery
acy(Bmat, imod.fast$M)
acy(Bmat, imod.imax$M)
acy(Bmat, imod.jade$M)

# compare source signal recovery
cor(Amat, imod.fast$S)
cor(Amat, imod.imax$S)
cor(Amat, imod.jade$S)


##########   EXAMPLE 2   ##########

# generate noiseless data (p != r)
set.seed(123)
nobs <- 1000
Amat <- cbind(icasamp("a", "rnd", nobs), icasamp("b", "rnd", nobs))
Bmat <- matrix(2 * runif(200), nrow = 100, ncol = 2)
Xmat <- tcrossprod(Amat, Bmat)

# ICA via different algorithms
imod.fast <- ica(Xmat, nc = 2)
imod.imax <- ica(Xmat, nc = 2, method = "imax")
imod.jade <- ica(Xmat, nc = 2, method = "jade")

# compare source signal recovery
cor(Amat, imod.fast$S)
cor(Amat, imod.imax$S)
cor(Amat, imod.jade$S)


##########   EXAMPLE 3   ##########

# generate noisy data (p != r)
set.seed(123)
nobs <- 1000
Amat <- cbind(icasamp("a", "rnd", nobs), icasamp("b", "rnd", nobs))
Bmat <- matrix(2 * runif(200), 100, 2)
Emat <- matrix(rnorm(10^5), nrow = 1000, ncol = 100)
Xmat <- tcrossprod(Amat,Bmat) + Emat

# ICA via different algorithms
imod.fast <- ica(Xmat, nc = 2)
imod.imax <- ica(Xmat, nc = 2, method = "imax")
imod.jade <- ica(Xmat, nc = 2, method = "jade")

# compare source signal recovery
cor(Amat, imod.fast$S)
cor(Amat, imod.imax$S)
cor(Amat, imod.jade$S)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
