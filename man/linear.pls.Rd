\name{linear.pls}
\Rdversion{1.1}
\alias{linear.pls}
\title{Linear Partial Least Squares
}
\description{This function computes the  Partial Least Squares solution and the first derivative of the regression coefficients.
}
\usage{
linear.pls(X, y, m = ncol(X),model.selection="aic")
}
\arguments{
 \item{X}{matrix of predictor observations.
}
  \item{y}{vector of response observations. The length of \code{y} is the same as the number of rows of \code{X}.
}
  \item{m}{maximal number of Partial Least Squares components. Default is \code{m}=ncol(X).
}
  \item{model.selection}{Which model selection criterion should be used? Element from \code{c("aic","bic","gmdl")}.
}
  
}
\details{We first standardize \code{X} to zero mean and unit variance.
}
\value{
The function returns an object of class "plsdof".
\item{coefficients}{matrix of regression coefficients}
\item{intercept}{vector of regression intercepts}
\item{DoF}{Degrees of Freedom}
\item{sigmahat}{vector of estimated model error}
\item{dBeta}{array of the first derivative of \code{coefficients}}
\item{covariance}{array of the covariance matrices of \code{coefficients}}
\item{m.opt}{optimal number of PLS components, as determined by \code{model.selection}}
}

\references{

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 

}
\author{Nicole Kraemer
}




\seealso{
\code{\link{kernel.pls.ic}}, \code{\link{kernel.pls.cv}},\code{\link{kernel.pls}}
}
\examples{
n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)


pls.object<-linear.pls(X,y,m=5,model.selection="bic")


}

\keyword{multivariate}
