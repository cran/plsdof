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
The function returns the variance-covariance matrix for the optimal number of components.
}
\value{
variance-covariance matrix
}
\references{
Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 


}
\author{
Nicole Kraemer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{coef.plsdof}}, \code{\link{linear.pls}}
}
\examples{

n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)


pls.object<-linear.pls(X,y,m=5,model.selection="bic")
myvcov<-vcov(pls.object)
}

\keyword{models}
