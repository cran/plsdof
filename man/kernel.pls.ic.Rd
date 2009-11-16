\name{kernel.pls.ic}
\Rdversion{1.1}
\alias{kernel.pls.ic}
\title{Model selection for Kernel Partial Least Squares based on information criteria
}
\description{This function computes the optimal model parameters using three different model selection criteria (aic, bic, gmdl) and based on
two different Degrees of Freedom estimates for PLS.
}
\usage{
kernel.pls.ic(X, y, m = ncol(X), type = "vanilla", sigma = 1, step.size = 1)
}
\arguments{
 \item{X}{matrix of predictor observations.
}
  \item{y}{vector of response observations. The length of \code{y} is the same as the number of rows of \code{X}.
}
  \item{m}{maximal number of Partial Least Squares components. Default is \code{m=ncol(X)}.
}
  \item{type}{type of kernel. \code{type="vanilla"} is a linear kernel. \code{type="gaussian"} is a gaussian kernel. Default is \code{type="vanilla"}.
}
  \item{sigma}{vector of kernel parameters. If \code{type="gaussian"}, these are the kernel widths. If the vanilla kernel is used,
  \code{sigma} is not used. Default value is \code{sigma=1}. 
}
  \item{step.size}{
After how many steps should the latent components be re-orthogonalized? See \code{\link{kernel.pls.fit}}  for more details. Default is \code{step.size=1}.
}
  
}
\details{For the linear kernel (\code{type="vanilla"}), we standardize \code{X} to zero mean and unit variance. For the Gaussian kernel (\code{type="gaussian"}), we normalize \code{X} such that the range of each column is [-1,1]. 

The default value for \code{sigma} is in general NOT a sensible parameter, and \code{sigma} should always be selected  from a RANGE of values. The default value for \code{m} is a sensible upper bound only for the vanilla kernel.
}
\value{
\item{DoF}{Degrees of Freedom}
\item{m.aic}{optimal number of components for aic}
\item{m.bic}{optimal number of components for bic}
\item{m.gmdl}{optimal number of components for gmdl}
\item{m.aic.naive}{optimal number of components for aic and the naive Degrees of Freedom}
\item{m.bic.naive}{optimal number of components for bic and the naive Degrees of Freedom}
\item{m.gmdl.naive}{optimal number of components for gmdl and the naive Degrees of Freedom}
\item{sigma.aic}{optimal sigma for aic, only returned if \code{type="gaussian"}}
\item{sigma.bic}{optimal sigma for bic, only returned if \code{type="gaussian"}}
\item{sigma.gmdl}{optimal sigma for gmdl, only returned if \code{type="gaussian"}}
\item{sigma.aic.naive}{optimal sigma for aic and the naive Degrees of Freedom, only returned if \code{type="gaussian"} }
\item{sigma.bic.naive}{optimal sigma for bic and the naive Degrees of Freedom, only returned if \code{type="gaussian"}}
\item{sigma.gmdl.naive}{optimal sigma for gmdl and the naive Degrees of Freedom, only returned if \code{type="gaussian"}}



}
\references{
Akaikie, H. (1973) "Information Theory and an Extension of the Maximum Likelihood Principle". Second International Symposium on Information Theory, 267 - 281.

Hansen, M., Yu, B. (2001). "Model Selection and Minimum Descripion Length Principle". Journal of
the American Statistical Association, 96, 746 - 774


Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 

Schwartz, G. (1979) "Estimating the Dimension of a Model" Annals of Statistics 26(5), 1651 - 1686.


}
\author{Nicole Kraemer, Mikio L. Braun
}

\seealso{
\code{\link{kernel.pls}}, \code{\link{kernel.pls.cv}}
}
\examples{
n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

# compute linear PLS
linear.pls<-kernel.pls.ic(X,y,m=ncol(X))
}

\keyword{multivariate}
