\name{kernel.pls}
\Rdversion{1.1}
\alias{kernel.pls}
\title{Kernel Partial Least Squares
}
\description{This function computes the Kernel Partial Least Squares fit.
}
\usage{
kernel.pls(X, y, m = ncol(X), Xtest, ytest, compute.DoF, type, sigma, step.size)
}
\arguments{
 \item{X}{matrix of predictor observations.
}
  \item{y}{vector of response observations. The length of \code{y} is the same as the number of rows of \code{X}.
}
  \item{m}{maximal number of Partial Least Squares components. Default is \code{m}=ncol(X).
}
  \item{Xtest}{optional matrix of test observations. Default is \code{Xtest=NULL}.
}
  \item{ytest}{optional vector of test observations. Default is \code{ytest=NULL}.
}
  \item{compute.DoF}{Logical variable. If \code{compute.DoF=TRUE}, the Degrees of Freedom of Partial Least Squares are computed. Default is \code{compute.DoF=FALSE}.
}
  \item{type}{type of kernel. \code{type="vanilla"} is a linear kernel. \code{type="gaussian"} is a gaussian kernel. Default is \code{type="vanilla"}.
}
  \item{sigma}{vector of kernel parameters. If \code{type="gaussian"}, these are the kernel widths. If the vanilla kernel is used,
  \code{sigma} is not used. Default value is \code{sigma}=1. 
}
  \item{step.size}{
After how many steps should the latent components be re-orthogonalized? See \code{\link{kernel.pls.fit}}  for more details. Default is \code{step.size=1}.
}
}
\details{For the linear kernel (\code{type="vanilla"}), we standardize \code{X} to zero mean and unit variance. For the Gaussian kernel (\code{type="gaussian"}), we normalize \code{X} such that the range of each column is [-1,1]. 

The default value for \code{sigma} is in general NOT a sensible parameter, and \code{sigma} should always be selected via a model selection criterion. The default value for \code{m} is a sensible upper bound only for the vanilla kernel.
}
\value{
\item{A}{array of kernel coefficients}
\item{prediction}{prediction on \code{Xtest}, if \code{Xtest} is provided}
\item{mse}{test error, if \code{Xtest} and \code{ytest} is provided}
\item{coefficients}{matrix of regression coefficients, if \code{type="vanilla"}}
\item{intercept}{vector of regression intercepts, if \code{type="vanilla"}}
\item{DoF}{Degrees of Freedom, if \code{compute.DoF=TRUE}}
\item{RSS}{residual sum of squares, if \code{compute.DoF=TRUE}}
\item{Yhat}{array of fitted values}
\item{sigmahat}{vector of estimated model error, if \code{compute.DoF=TRUE}}
\item{yhat}{sum of squared fitted values}
}

\references{

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 

}
\author{Nicole Kraemer, Mikio L. Braun
}




\seealso{
\code{\link{kernel.pls.ic}}, \code{\link{kernel.pls.cv}}
}
\examples{
n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

ntest<-200 #
Xtest<-matrix(rnorm(ntest*p),ncol=p) # test data
ytest<-rnorm(ntest) # test data

# compute linear PLS + degrees of freedom + prediction on Xtest
linear.object<-kernel.pls(X,y,m=ncol(X),compute.DoF=TRUE,Xtest=Xtest,ytest=NULL)

# compute nonlinear PLS + test error
sigma=exp(seq(0,4,length=10))
nonlinear.object=kernel.pls(X,y,m=10,type="gaussian",sigma=sigma,Xtest=Xtest,ytest=ytest)
}

\keyword{multivariate}
