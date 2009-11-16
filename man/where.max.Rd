\name{where.max}
\Rdversion{1.1}
\alias{where.max}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute argmax of a vector
}
\description{This function computes the argmax of a matrix.
}
\usage{
where.max(M)
}
\arguments{
  \item{M}{matrix
}
}
\details{The function returns the row and the column index of the argmax of M.
}
\value{vector that contains the row index and the column index
}

\author{Nicole Kraemer
}


\examples{
M<-matrix(sample(20),ncol=4)
index<-which.max(M)
}
\keyword{math}