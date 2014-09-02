\name{icasamp}
\alias{icasamp}
\title{
Sample from various source signal distributions
}
\description{
Sample observations from the 18 source signal distributions used in Bach and Jordan (2002). Can also return density values and kurtosis for each distribution. Use \code{\link{icaplot}} to plot distributions.
}
\usage{
icasamp(dname,query=c("rnd","pdf","kur"),nsamp=NULL,data=NULL)
}
\arguments{
  \item{dname}{
Distribution name: letter "a" through "r" (see Bach & Jordan, 2002).
}
  \item{query}{
What to return: \code{query="rnd"} for random sample, \code{query="pdf"} for density values, and \code{query="kur"} for kurtosis.
}
  \item{nsamp}{
Number of observations to sample. Only used if \code{query="rnd"}.
}
  \item{data}{
Data values for density evaluation. Only used if \code{query="pdf"}.
}
}
\value{
If \code{query="rnd"}, returns random sample of size \code{nsamp}.

If \code{query="pdf"}, returns density for input \code{data}.

If \code{query="kur"}, returns kurtosis of distribution.
}
\references{
Bach, F.R. (2002). \emph{kernel-ica}. MATLAB toolbox (ver 1.2) http://www.di.ens.fr/~fbach/kernel-ica/.
  
Bach, F.R. & Jordan, M.I. (2002). Kernel independent component analysis. \emph{Journal of Machine Learning Research, 3}, 1-48.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
Inspired by \code{usr_distrib.m} from Bach's (2002) \code{kernel-ica} MATLAB toolbox.
}
\examples{
##########   EXAMPLE   ##########

# sample 1000 observations from distribution "f"
set.seed(123)
mysamp<-icasamp("f","rnd",nsamp=1000)
xr<-range(mysamp)
hist(mysamp,freq=FALSE,ylim=c(0,.8),breaks=sqrt(1000))

# evaluate density of distribution "f"
xseq<-seq(-5,5,length.out=1000)
mypdf<-icasamp("f","pdf",data=xseq)
lines(xseq,mypdf)

# evaluate kurtosis of distribution "f"
icasamp("f","kur")

}
