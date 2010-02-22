\name{pls.cv}
\Rdversion{1.1}
\alias{pls.cv}
\title{Model selection for Partial Least Squares based on cross-validation
}
\description{This function computes the optimal model parameter using cross-validation.
}
\usage{
pls.cv(X, y, k, m,use.kernel=FALSE,compute.covariance=FALSE)
}
\arguments{
 \item{X}{matrix of predictor observations.
}
  \item{y}{vector of response observations. The length of \code{y} is the same as the number of rows of \code{X}.
}
\item{k}{number of cross-validation splits. Default is 10.}
  \item{m}{maximal number of Partial Least Squares components. Default is \code{m=ncol(X)}.
}
  \item{use.kernel}{Use kernel representation? Default is \code{use.kernel=FALSE}.
}
  
\item{compute.covariance}{If \code{TRUE}, the function computes the covariance for the cv-optimal regression coefficients.} 
  
}
\details{to do
}
\value{
The function returns an object of class "plsdof".
\item{cv.error}{vector of cross-validated errors}
\item{m.opt}{optimal number of components}
\item{intercept}{intercept}
\item{coefficients}{vector of regression coefficients}
\item{covariance}{If \code{TRUE} and \code{use.kernel=FALSE}, the covariance of the cv-optimal regression coefficients is returned.}


}
\references{

Kraemer, N., Sugiyama M. (2010). "The Degrees of Freedom of Partial Least Squares Regression". preprint, \url{http://arxiv.org/abs/1002.4112}



Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 


}
\author{Nicole Kraemer, Mikio L. Braun
}

\seealso{
\code{\link{pls.model}}, \code{\link{pls.ic}}
}
\examples{
n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

# compute linear PLS
pls.object<-pls.cv(X,y,m=ncol(X))
}

\keyword{multivariate}
