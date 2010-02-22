\name{benchmark.pls}
\Rdversion{1.1}
\alias{benchmark.pls}
\title{Comparison of model selection criteria for Partial Least Squares Regression.
}
\description{This function computes the test error over several runs  for different model selection strategies.
}
\usage{
benchmark.pls(X,y,m,R,ratio,verbose,k,ratio.samples,use.kernel=FALSE,criterion="bic")
}
\arguments{
 \item{X}{matrix of predictor observations.
}
  \item{y}{vector of response observations. The length of \code{y} is the same as the number of rows of \code{X}.
}
  \item{m}{maximal number of Partial Least Squares components. Default is \code{m=ncol(X)}.
}
\item{R}{number of runs. Default is 20.}
\item{ratio}{ratio  no of training examples/(no of training examples + no of test examples). Default is 0.8}
\item{verbose}{If \code{TRUE}, the functions plots the progress of the function. Default is \code{TRUE}.}
\item{k}{number of cross-validation splits. Default is 10.}
\item{ratio.samples}{Ratio of (no of training examples + no of test examples)/\code{nrow(X)}. Default is 1. }
\item{use.kernel}{Use kernel representation? Default is \code{use.kernel=FALSE}.}
  \item{criterion}{Choice of the model selection criterion. One of the three options aic, bic, gmdl. Default is "bic".
}
  

  
}
\details{The function estimates the optimal number of PLS components based on four different criteria: (1) Cross-validation,
(2) information criteria with the naive Degrees of Freedom DoF(m)=M+1, (3) information criteri with the Degrees of Freedom computed via
a Lanczos represenation of PLS and (4)  information criteri with the Degrees of Freedom computed via
a Krylov represenation of PLS. Note that the latter two options only differ with respect to the estimation of the model error.

In addition, the function computes the test error of the "zero model", i.e. \code{mean(y)} on the training data is used for prediction.
}
\value{
\item{MSE}{data frame of size R x 5. It contains the test error for the five different methods over the R runs.}
\item{M}{data frame of size R x 5. It contains the optimal number of components for the five different methods over the R runs.}
\item{DoF}{data frame of size R x 5. It contains the optimal Degrees of Freedom for the five different methods over the R runs.}
\item{TIME}{data frame of size R x 4. It contains the run for all methods (apart from the zero model) over the R runs.}
\item{M.CRASH}{data frame of size R x 2. It contains the number of components for which the Krylov representation and the Lanczos representation return negative Degrees of Freedom, hereby indicating numerical problems.}

}
\references{

Kraemer, N., Sugiyama M. (2010). "The Degrees of Freedom of Partial Least Squares Regression". preprint, \url{http://arxiv.org/abs/1002.4112}




}
\author{Nicole Kraemer
}

\seealso{
\code{\link{pls.ic}}, \code{\link{pls.cv}}
}
\examples{
# to do
}

\keyword{multivariate}
