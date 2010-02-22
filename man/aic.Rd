\name{aic}
\Rdversion{1.1}
\alias{aic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Akaike Information Criterion
}
\description{
This function computes the Akaike Information Criterion and its argmin.}
\usage{
aic(RSS, n, DoF, sigmahat)
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
The aic criterion is defined as 
\deqn{{bic}= \frac{{RSS}}{n} + 2\frac{{DoF}}{n} \sigma^ 2\,.}
Note that it is also possible to use this function for other regression methods than Partial Least Squares.
}

\value{
\item{score}{the vector of the aic values}
\item{par}{index of the first local minimum of \code{score}}
}

\references{
Akaikie, H. (1973) "Information Theory and an Extension of the Maximum Likelihood Principle".
Second International Symposium on Information Theory, 267 - 281.

Kraemer, N., Sugiyama M. (2010). "The Degrees of Freedom of Partial Least Squares Regression". preprint, \url{http://arxiv.org/abs/1002.4112}

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 

 

}
\author{Nicole Kraemer, Mikio L. Braun
}
\seealso{
 \code{\link{information.criteria}}, \code{\link{pls.ic}}, \code{\link{bic}}, \code{\link{gmdl}}
}
\examples{
## this is an internal function called by information.criteria() and pls.ic()
}

\keyword{models}
