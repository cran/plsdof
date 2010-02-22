\name{vcov.plsdof}
\Rdversion{1.1}
\alias{vcov.plsdof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variance-covariance matrix
}
\description{
This function returns the variance-covariance matrix of a plsdof-object.
}
\usage{
\method{vcov}{plsdof}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class "plsdof" that is returned by the function \code{linear.pls}
}
\item{...}{additional parameters}
}
\details{
The function returns the variance-covariance matrix for the optimal number of components. It can be applied to objects returned by
\code{pls.ic} and \code{pls.cv}.
}
\value{
variance-covariance matrix
}
\references{
Kraemer, N., Sugiyama M. (2010). "The Degrees of Freedom of Partial Least Squares Regression". preprint, \url{http://arxiv.org/abs/1002.4112}



Kraemer, N., Sugiyama M., Braun, M.L. (2009) "Lanczos Approximations for the Speedup of Kernel Partial Least Squares Regression." Proceedings of the Twelfth International Conference on Artificial Intelligence and Statistics (AISTATS), p. 272-279

}
\author{
Nicole Kraemer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{coef.plsdof}}, \code{\link{pls.ic}}, \code{\link{pls.cv}}
}
\examples{

n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)


pls.object<-pls.ic(X,y,m=5,criterion="bic")
myvcov<-vcov(pls.object)
}

\keyword{models}
