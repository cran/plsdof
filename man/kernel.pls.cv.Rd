\name{kernel.pls.cv}
\Rdversion{1.1}
\alias{kernel.pls.cv}
\title{Cross-validation for Kernel Partial Least Squares
}
\description{
This function computes the cross-validation-optimal model parameters for Kernel Partial Least Squares.
}
\usage{
kernel.pls.cv(X, y, k = 10, m = ncol(X), type = "vanilla", sigma = 1)
}
\arguments{
  \item{X}{matrix of predictor observations.
}
  \item{y}{vector of response observations. The length of \code{y} is the same as the number of rows of \code{X}.
}
  \item{k}{number of cross-validation splits. Default is \code{k}=10.}
  \item{m}{maximal number of Partial Least Squares components. Default is \code{m}=ncol(X).}
  Default is \code{type="vanilla"}.
\item{type}{type of kernel. \code{type="vanilla"} is a linear kernel. \code{type="gaussian"} is a gaussian kernel. Default is \code{type="vanilla"}.
}
  \item{sigma}{vector of kernel parameters. If \code{type="gaussian"}, these are the kernel widths. If the vanilla kernel is used,
  \code{sigma} is not used. Default value is \code{sigma}=1. 
}
}
\details{For the linear kernel (\code{type="vanilla"}), we standardize \code{X} to zero mean and unit variance. For the Gaussian kernel (\code{type="gaussian"}), we normalize \code{X} such that the range of each column is [-1,1]. 

The default value for \code{sigma} is in general NOT a sensible parameter, and \code{sigma} should always be selected via cross-validation from a RANGE of values. The default value for \code{m} is a sensible upper bound only for the vanilla kernel.

}
\value{
\item{cv.error}{cross-validated error. If \code{type="vanilla"}, this is a vector of length \code{m}. If \code{type="gaussian"},
this is a matrix of dimensions \code{length(sigma)} x \code{m}.}
\item{m.opt}{optimal numnber of Partial Least Squares components}
\item{sigma.opt}{optimal kernel paramter. This is only returned if \code{type="gaussian"}.}

}

\author{Nicole Kraemer, Mikio L. Braun
}


\seealso{
\code{\link{kernel.pls.ic}}, \code{\link{kernel.pls}}
}
\examples{
n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

# compute linear PLS
linear.pls<-kernel.pls.cv(X,y,m=ncol(X),k=5)

# compute nonlinear PLS
sigma=exp(seq(0,4,length=10))
nonlinear.pls=kernel.pls(X,y,m=10,type="gaussian",sigma=sigma)
}

\keyword{multivariate}
