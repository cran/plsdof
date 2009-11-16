\name{information.criteria}
\Rdversion{1.1}
\alias{information.criteria}
\title{Information criteria
}
\description{
This function computes the optimal model parameters using three different model selection criteria (aic, bic, gmdl) and based on
two different Degrees of Freedom estimates for Partial Least Squares.
}
\usage{
information.criteria(RSS, DoF, yhat,sigmahat, n, DoF.max)
}
\arguments{
  \item{RSS}{matrix of residual sum of squares. The columns correspond to the number of latent components for PLS, and
  the rows correspond to the different kernel parameters. If the vanilla kernel is used, RSS consists of one row. For a Gaussian
  kernel, each row corresponds to a sigma value.
}
  \item{DoF}{
Degrees of Freedom. The size of \code{DoF} is the same as the size of \code{RSS}.
}
  \item{yhat}{matrix of squared norm of yhat. The size of \code{yhat} is the same as the size of \code{RSS}.
}
  \item{sigmahat}{Estimated model error. The size of \code{sigmahat} is the same as the size of \code{RSS}.
}
  \item{n}{number of observations.
}
  \item{DoF.max}{Upper bound for the degrees of freedom.
}
}
\details{In an initial step, we first normalize the \code{DoF} to lie in the range 0 and \code{DoF.max}. For linear Partial Least Squares, a sensible upper bound \code{DoF.max} is given by \code{DoF.max=min(nrow(X)-1,ncol(X))}.}
\value{
 \item{DoF}{degrees of freedom}
\item{par.aic}{indices for the aic-optimal model}
\item{par.bic}{indices for the bic-optimal model}
\item{par.gmdl}{indices for the gmdl-optimal model}
\item{par.aic.naive}{indices for the aic-optimal model if the naive Degrees of Freedom are used}
\item{par.bic.naive}{indices for the bic-optimal model if the naive Degrees of Freedom are used}
\item{par.gmdl.naive}{indices for the gmdl-optimal model if the naive Degrees of Freedom are used}

}
\references{
Akaikie, H. (1973) "Information Theory and an Extension of the Maximum Likelihood Principle". Second International Symposium on Information Theory, 267 - 281.

Hansen, M., Yu, B. (2001). "Model Selection and Minimum Descripion Length Principle". Journal of
the American Statistical Association, 96, 746 - 774


Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 

Schwartz, G. (1979) "Estimating the Dimension of a Model" Annals of Statistics 26(5), 1651 - 1686.

}
\author{
Nicole Kraemer, Mikio Braun
}

\seealso{
\code{\link{kernel.pls.ic}}
}
\examples{
## This is an internal function called by kernel.pls.ic

}

\keyword{model}
