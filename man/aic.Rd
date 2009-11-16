\name{aic}
\Rdversion{1.1}
\alias{aic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Akaike information criterion
}
\description{
This function computes the optimal parameter set based on the Akaike information criterion.
}
\usage{
aic(RSS, n, DoF, sigmahat)
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
The aic criterion is defined as 
\deqn{{aic}= \frac{{RSS}}{n} + 2\frac{{DoF}}{n} \sigma^ 2\,.}
For the aic criterion, we need an estimtate for the model error \eqn{\sigma}. We use \code{min(sigmahat)} as an estimate. Note that it is also possible to use this function for other regression methods than Partial Least Squares.
}
\value{
The function returns a vector of length 2 that contains the row and column index of the aic-optimal model.
}
\references{
Akaikie, H. (1973) "Information Theory and an Extension of the Maximum Likelihood Principle". Second International Symposium on Information Theory, 267 - 281.

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 


}
\author{
Nicole Kraemer, Mikio L. Braun
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{information.criteria}}, \code{\link{kernel.pls.ic}}, \code{\link{bic}}, \code{\link{gmdl}}
}
\examples{
## this is an internal function called by information.criteria() and kernel.pls.ic()
}

\keyword{models}
