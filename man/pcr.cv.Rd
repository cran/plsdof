\name{pcr.cv}
\Rdversion{1.1}
\alias{pcr.cv}
\title{Model selection for Princinpal Components regression based on cross-validation
}
\description{This function computes the optimal model parameter using cross-validation.
}
\usage{
pcr.cv(X,y,k=10,groups=NULL,scale=TRUE,eps=0.000001,plot.it=FALSE)
}
\arguments{
 \item{X}{matrix of predictor observations.
}
  \item{y}{vector of response observations. The length of \code{y} is the same as the number of rows of \code{X}.
}
\item{k}{number of cross-validation splits. Default is 10.}
\item{groups}{an optional vector with the same length as \code{y}. It encodes a partitioning of the data into distinct subgroups. If \code{groups} is provided, \code{k=10} is ignored and instead, cross-validation is performed based on the partioning. Default is \code{NULL}.}
\item{scale}{Should the predictor variables be scaled to unit variance? Default is \code{TRUE}.}
  \item{eps}{precision. Eigenvalues of the correlation matrix of \code{X} that are smaller than \code{eps} are set to 0. The default value is \code{eps=10^{-6}.}}
  \item{plot.it}{Logical. If \code{TRUE}, the function plots the cross-validation-error as a function of the number of components. Default is \code{FALSE}.}
}
\details{The function computes the principal components on the scaled predictors. Based on the regression coefficients \code{coefficients.jackknife} computed on the cross-validation
splits, we can estimate their mean and their variance using the jackknife. We remark that under a fixed design and the assumption of
normally distributed \code{y}-values, we can also derive the true distribution of the regression coefficients.
}
\value{
The function returns an object of class "plsdof".
\item{cv.error.matrix}{matrix of cross-validated errors. A row corresponds to one cross-validation split.}
\item{cv.error}{vector of cross-validated errors}
\item{m.opt}{optimal number of components}
\item{intercept}{intercept}
\item{coefficients}{vector of regression coefficients}
\item{coefficients.jackknife}{Array of the regression coefficients on each of the cross-validation splits. The dimension is \code{ncol(X) x (ncol(X)+1) x k}.}


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

# compute PCR 
pcr.object<-pcr.cv(X,y,scale=FALSE)
pcr.object1<-pcr.cv(X,y,groups=sample(c(1,2,3),n,replace=TRUE))
}

\keyword{multivariate}
