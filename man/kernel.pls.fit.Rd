\name{kernel.pls.fit}
\Rdversion{1.1}
\alias{kernel.pls.fit}
\title{Kernel Partial Least Squares fit
}
\description{
This function computes the Kernel Partial Least Squares fit.
}
\usage{
kernel.pls.fit(K, y, m, compute.DoF, step.size)
}
\arguments{
  \item{K}{kernel matrix. This is a square matrix with dimension = number of observations. We assume that the data is centered
  in the feature space.
}
  \item{y}{vector of response observations. 
}
  \item{m}{number of Partial Least Squares components
}
  \item{compute.DoF}{Logical variable. If \code{compute.DoF=TRUE}, the Degrees of Freedom of Partial Least Squares are computed.
}
  \item{step.size}{After how many steps should the latent components be re-orthogonalized? See below for more details.
}
}
\details{This function computed the kernel coefficients for Partial Least Squares and its Degrees of Freedom. The computation
relies on a sparse matrix decomposition of the predictor variables. The sparse structure can sometimes lead to non-orthogonal latent components
(due to numerical problems). To avoid this, we orthogonalize the latent components after \code{step.size} iterations.
}
\value{
\item{Alpha}{matrix of kernel coefficients. This is a matrix of size \code{length(y)} x \code{m}.}
\item{dYhat}{array of first derivatives of \eqn{\widehat y} with respect to \code{y}. The dimension of \code{dYhat} equals
[m,length(y),length(y)].}
\item{Yhat}{matrix of fitted values. This is a matrix of size \code{length(y)} x \code{m}.}
\item{DoF}{vector of Degrees of Freedom.}
\item{RSS}{vector of residual sum of squares}
\item{sigmahat}{vector of estimated model error}
}

\references{

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 

}
\author{Nicole Kraemer, Mikio L. Braun
}

\seealso{
\code{\link{kernel.pls}}, \code{\link{kernel.pls.cv}}, \code{\link{kernel.pls.ic}}
}
\examples{
## this is an internal function called by kernel.pls
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
