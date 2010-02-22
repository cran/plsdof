\name{first.local.minimum}
\Rdversion{1.1}
\alias{first.local.minimum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Index of the first local minimum.
}
\description{
This function computes the index of the first local minimum.
}
\usage{
first.local.minimum(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector.

}
  
}

\value{
the index of the first local minimum of \code{x}.
}


\references{
Kraemer, N., Sugiyama M. (2010). "The Degrees of Freedom of Partial Least Squares Regression". preprint, \url{http://arxiv.org/abs/1002.4112}


}
\author{
Nicole Kraemer
}


\examples{
v<-rnorm(30)
out<-first.local.minimum(v)
}
\keyword{math}
