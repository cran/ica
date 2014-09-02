\name{icaplot}
\alias{icaplot}
\title{
Plot densities of source signal distributions
}
\description{
Plot density (pdf) and kurtosis for the 18 source signal distributions used in Bach and Jordan (2002); see \code{\link{icasamp}} for more information.
}
\usage{
icaplot(xseq=seq(-2,2,length.out=500),xlab="",ylab="",
        lty=1,lwd=1,col="black",...)
}
\arguments{
  \item{xseq}{
Sequence of ordered data values for plotting density.
}
  \item{xlab}{
X-axis label for plot (default is no label).
}
  \item{ylab}{
Y-axis label for plot (default is no label).
}
  \item{lty}{
Line type for each density (scalar or vector of length 18).
}
  \item{lwd}{
Line width for each density (scalar or vector of length 18).
}
  \item{col}{
Line color for each density (scalar or vector of length 18).
}
  \item{...}{
Optional inputs for \code{plot}.
}
}
\value{
Produces a plot with \code{NULL} return value.
}
\references{
Bach, F.R. (2002). \emph{kernel-ica}. MATLAB toolbox (ver 1.2) http://www.di.ens.fr/~fbach/kernel-ica/.
  
Bach, F.R. & Jordan, M.I. (2002). Kernel independent component analysis. \emph{Journal of Machine Learning Research, 3}, 1-48.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\examples{
##########   EXAMPLE (not run)   ##########

#quartz(height=9,width=7)
#par(mar=c(3,3,3,3))
#icaplot()

}
