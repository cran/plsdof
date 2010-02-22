\name{gmdl}
\Rdversion{1.1}
\alias{gmdl}
\title{Generalized minimum description length
}
\description{
This function computes the generalized minimum description length and its argmin.
}
\usage{
gmdl(sigmahat, n, DoF, yhat)
}
\arguments{
   \item{sigmahat}{vector of estimated model errors.
}
  \item{n}{number of observations.
}
  \item{DoF}{vector of Degrees of Freedom. The length of \code{DoF} is the same as the length of \code{sigmahat}.
}
  \item{yhat}{vector of squared norm of yhat. The length of \code{yhat} is the same as the length of \code{sigmahat}
}
}

\details{
The gmdl criterion is defined as
\deqn{gmdl=\frac{n}{2}log(S)+\frac{DoF}{2}log(F)+\frac{1}{2}log(n)}
with \deqn{S=sigmahat^2}
Note that it is also possible to use this function for other regression methods than Partial Least Squares.


}

\value{
\item{score}{the vector of the aic values}
\item{par}{index of the first local minimum of \code{score}}
}
\references{
Hansen, M., Yu, B. (2001). "Model Selection and Minimum Descripion Length Principle". Journal of
the American Statistical Association, 96, 746 - 774

Kraemer, N., Sugiyama M. (2010). "The Degrees of Freedom of Partial Least Squares Regression". preprint, \url{http://arxiv.org/abs/1002.4112}


Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 

}
\author{Nicole Kraemer, Mikio L. Braun
}


\seealso{
\code{\link{information.criteria}}, \code{\link{pls.ic}}, \code{\link{aic}}, \code{\link{bic}}
}
\examples{
## This is an internal function.
}
\keyword{math}
