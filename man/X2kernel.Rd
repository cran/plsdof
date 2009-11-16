\name{X2kernel}
\Rdversion{1.1}
\alias{X2kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of the kernel matrix}
\description{This function computes the Gram matrices for a gaussian kernel and differnt kernel widths.
}
\usage{
X2kernel(X, Xtest = NULL, sigma = 1)
}
\arguments{
 \item{X}{matrix of predictor observations.
}

  \item{Xtest}{optional matrix of test observations. Default is \code{Xtest=NULL}.
}
  \item{sigma}{vector of kernel widths. Default is \code{sigma=1}.
}
}
\details{We first scale \code{X} such that the range of each column lies in [-1,1]. Additionally, the data is centered in feature space.
If test data is provided, the same transformations are also applied to them.

The default value for \code{sigma} is in general NOT a sensible parameter, and \code{sigma} should always be selected via a model selection criterion.
}
\value{
    \item{K}{array of kernel matrices}
    \item{Ktest}{array of kernel matrices for \code{Xtest}, if provided}
    \item{sigma}{vector of kernel widths}
}

\author{Nicole Kraemer
}

\seealso{
\code{\link{X2kernel}}, \code{\link{kernel.pls}}, \code{\link{kernel.pls.ic}}, \code{\link{kernel.pls.cv}}
}
\examples{
p<-20 # number of variables
n<-100 # number of observations
ntest<-50 # number of test observations

X<-matrix(rnorm(n*p),ncol=p)
Xtest<-matrix(rnorm(ntest*p),ncol=p)
sigma<-exp(seq(0,4,length=10))

kernel.object<-X2kernel(X,Xtest,sigma=sigma)

}

\keyword{multivariate}
