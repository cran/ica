\name{ica-package}
\alias{ica-package}
\alias{ica}
\docType{package}
\title{Independent Component Analysis}
\description{
Independent Component Analysis (ICA) using various algorithms: FastICA, Information-Maximization (Infomax), and Joint Approximate Diagonalization of Eigenmatrices (JADE).
}
\details{
%\tabular{ll}{
%Package: \tab ica\cr
%Type: \tab Package\cr
%Version: \tab 1.0\cr
%Date: \tab 2013-08-31\cr
%License: \tab GPL (>=2)\cr
%}
The functions \code{\link{icafast}}, \code{\link{icaimax}}, and \code{\link{icajade}} calculate ICA demcompositions using the FastICA, Infomax, and JADE algorithms (respectively). The function \code{\link{icasamp}} can be used to sample from various interesting distirubtions, which are useful for comparing ICA algorithms.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>

Maintainer: Nathaniel E. Helwig <helwig@umn.edu>
}
\references{
  Amari, S., Cichocki, A., & Yang, H.H. (1996). A new learning algorithm for blind signal separation. In D. S. Touretzky, M. C. Mozer, and M. E. Hasselmo (Eds.), \emph{Advances in Neural Information Processing Systems, 8}. Cambridge, MA: MIT Press.

  Bach, F.R. (2002). \emph{kernel-ica}. MATLAB toolbox (ver 1.2) http://www.di.ens.fr/~fbach/kernel-ica/.
  
  Bach, F.R. & Jordan, M.I. (2002). Kernel independent component analysis. \emph{Journal of Machine Learning Research, 3}, 1-48.
  
  Bell, A.J. & Sejnowski, T.J. (1995). An information-maximization approach to blind separation and blind deconvolution. \emph{Neural Computation, 7}, 1129-1159.
  
  Cardoso, J.F., & Souloumiac, A. (1993). Blind beamforming for non-Gaussian signals. \emph{IEE Proceedings-F, 140}, 362-370.
  
  Cardoso, J.F., & Souloumiac, A. (1996). Jacobi angles for simultaneous diagonalization. \emph{SIAM Journal on Matrix Analysis and Applications, 17}, 161-164.
  
  Helwig, N.E. (in prep). On the relationship between FastICA and Infomax: Fast and robust fixed point algorithms for information-maximization.
  
  Helwig, N.E. & Hong, S. (2013). A critique of Tensor Probabilistic Independent Component Analysis: Implications and recommendations for multi-subject fMRI data analysis. \emph{Journal of Neuroscience Methods, 213}, 263-273.
  
  Hyvarinen, A. (1999). Fast and robust fixed-point algorithms for independent component analysis. \emph{IEEE Transactions on Neural Networks, 10}, 626-634.
  
  Tucker, L.R. (1951). \emph{A method for synthesis of factor analysis studies} (Personnel Research Section Report No. 984). Washington, DC: Department of the Army.
}
\keyword{ package }
\examples{
# See examples for icafast, icaimax, icajade, and icasamp
}
