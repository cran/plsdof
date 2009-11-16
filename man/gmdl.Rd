\name{gmdl}
\Rdversion{1.1}
\alias{gmdl}
\title{Generalized minimum description length
}
\description{
This function computes the optimal parameter set based on the generalized minimum description length.
}
\usage{
gmdl(RSS, n, DoF, yhat)
}
\arguments{
  \item{RSS}{matrix of residual sum of squares. The columns correspond to the number of latent components for PLS, and
  the rows correspond to the different kernel parameters. If the vanilla kernel is used, RSS consists of one row. For a Gaussian
  kernel, each row corresponds to a sigma value.
}
  \item{n}{number of observations.
}
  \item{DoF}{Degrees of Freedom. The size of \code{DoF} is the same as the size of \code{RSS}.
}
  \item{yhat}{matrix of squared norm of yhat. The size of \code{yhat} is the same as the size of \code{RSS}
}
}

\details{
The gmdl criterion is defined as
\deqn{gmdl=\frac{n}{2}log(S)+\frac{DoF}{2}log(F)+\frac{1}{2}log(n)}
with \deqn{S=\frac{RSS}{n-DoF}} and\deqn{F=\frac{\|\hat y\|^ 2}{DoF\cdot S}}
Note that it is also possible to use this function for other regression methods than Partial Least Squares.


}

\value{
The function returns a vector of length 2 that contains the row and column index of the gmdl-optimal model
}
\references{
Hansen, M., Yu, B. (2001). "Model Selection and Minimum Descripion Length Principle". Journal of
the American Statistical Association, 96, 746 - 774

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 

}
\author{Nicole Kraemer, Mikio L. Braun
}


\seealso{
\code{\link{information.criteria}}, \code{\link{kernel.pls.ic}}, \code{\link{aic}}, \code{\link{bic}}
}
\examples{
## This is an internal function.
}
\keyword{math}
