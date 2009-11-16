\name{bic}
\Rdversion{1.1}
\alias{bic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian information criterion
}
\description{
This function computes the optimal parameter set based on the Bayesian information criterion.
}
\usage{
bic(RSS, n, DoF, sigmahat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RSS}{matrix of residual sum of squares. The columns correspond to the number of latent components for PLS, and
  the rows correspond to the different kernel parameters. If the vanilla kernel is used, RSS consists of one row. For a Gaussian
  kernel, each row corresponds to a kernel width.
}
  \item{n}{number of observations.
}
  \item{DoF}{Degrees of Freedom. The size of \code{DoF} is the same as the size of \code{RSS}.
}
  \item{sigmahat}{Estimated model error. The size of \code{sigmahat} is the same as the size of \code{RSS}.
}
}
\details{
The bic criterion is defined as 
\deqn{{bic}= \frac{{RSS}}{n} + log(n)\frac{{DoF}}{n} \sigma^ 2\,.}
For the bic criterion, we need an estimtate for the model error \eqn{\sigma}. We use \code{min(sigmahat)} as an estimate. Note that it is also possible to use this function for other regression methods than Partial Least Squares.
}
\value{
The function returns a vector of length 2 that contains the row and column index of the bic-optimal model.
}
\references{
Schwartz, G. (1979) "Estimating the Dimension of a Model" Annals of Statistics 26(5), 1651 - 1686.

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 

 

}
\author{Nicole Kraemer, Mikio L. Braun
}
\seealso{
 \code{\link{information.criteria}}, \code{\link{kernel.pls.ic}}, \code{\link{aic}}, \code{\link{gmdl}}
}
\examples{
## this is an internal function called by information.criteria() and kernel.pls.ic()
}

\keyword{models}
