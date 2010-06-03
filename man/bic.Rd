\name{bic}
\Rdversion{1.1}
\alias{bic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian information criterion
}
\description{
This function computes the Bayesian Information Criterion and its argmin.
}
\usage{
bic(RSS, n, DoF, sigmahat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RSS}{vector of residual sum of squares.
}
  \item{n}{number of observations.
}
  \item{DoF}{vector of Degrees of Freedom. The length of \code{DoF} is the same as the length of \code{RSS}.
}
  \item{sigmahat}{Estimated model error. The length of \code{sigmahat} is the same as the length of \code{RSS}.
}
}
\details{
The bic criterion is defined as 
\deqn{{bic}= \frac{{RSS}}{n} + log(n)\frac{{DoF}}{n} \sigma^ 2\,.}
Note that it is also possible to use this function for other regression methods than Partial Least Squares.
}

\value{
\item{score}{the vector of the bic values}
\item{par}{index of the first local minimum of \code{score}}
}
\references{
Schwartz, G. (1979) "Estimating the Dimension of a Model" Annals of Statistics 26(5), 1651 - 1686.

Kraemer, N., Sugiyama M. (2010). "The Degrees of Freedom of Partial Least Squares Regression". preprint \url{http://arxiv.org/abs/1002.4112}


Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 

 

}
\author{Nicole Kraemer, Mikio L. Braun
}
\seealso{
 \code{\link{information.criteria}}, \code{\link{pls.ic}}, \code{\link{aic}}, \code{\link{gmdl}}
}
\examples{
## this is an internal function called by information.criteria() and kernel.pls.ic()
}

\keyword{models}
