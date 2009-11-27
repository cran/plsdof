\name{coef.plsdof}
\Rdversion{1.1}
\alias{coef.plsdof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Regression coefficients
}
\description{
This function returns the regression coefficients of a plsdof-object.
}
\usage{
\method{coef}{plsdof}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class "plsdof" that is returned by the function \code{linear.pls}
}
\item{...}{additional parameters}
}
\details{
The function returns the regression coefficients (without intercept) for the optimal number of components.
}
\value{
regression coefficients.
}
\references{

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 


}
\author{
Nicole Kraemer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{vcov.plsdof}}, \code{\link{linear.pls}}
}
\examples{

n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)


pls.object<-linear.pls(X,y,m=5,model.selection="bic")
mycoef<-coef(pls.object)
}

\keyword{models}
