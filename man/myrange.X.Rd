\name{myrange.X}
\Rdversion{1.1}
\alias{myrange.X}
\title{
Scaling of predictor variables
}
\description{This function scales each column of the matrix of predictor observations via a linear transformation ax+b.
}
\usage{
myrange.X(X, a = NULL, b = NULL)
}
\arguments{
  \item{X}{matrix of predictor obersvations
}
  \item{a}{vector of slopes for the linear transformation 
}
  \item{b}{vector of intercepts for the linear transformation
}
}
\details{If both \code{a} and \code{b} are given, the function transforms the ith column of \code{X} via \code{a[i]*X[,i]+b[i]}. 
Otherwise, the columns are scaled such that the range of each column is [-1,+1]. 
}
\value{
\item{X}{scaled X}
\item{a}{vector of slopes}
\item{b}{vector of intercepts}
}

\author{Nicole Kraemer, Mikio L. Braun
}

\seealso{
\code{\link{X2kernel}}
}
\examples{
p<-20 # number of variables
n<-100 # number of observations
ntest<-50 # number of test observations

X<-matrix(rnorm(n*p),ncol=p)
Xtest<-matrix(rnorm(ntest*p),ncol=p)

# scale X
range.object<-myrange.X(X)
X.scale<-range.object$X
# scale test data
Xtest.scale<-myrange.X(Xtest,a=range.object$a,b=range.object$b)$X

}

\keyword{math}
