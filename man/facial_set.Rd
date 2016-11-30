\name{facial_set}
\alias{facial_set}
\Rdversion{1.0.1}
\title{Finds the facial set}
\description{Identifies sampling zeros that contribute to the non-existence of the MLE. These problematic zero cells can be removed from the contingency table and the model can then be fit (as far as is possible) using the glm package.}
\usage{facial_set (data, formula)}
\alias{eMLEloglin}
\arguments{
  \item{data}{A data frame containing the contingency table. All cells must be included in data and the last column must be the cell counts. The number of variables in the contingency table must be between 2 and 8.}
  \item{formula}{A model formula.}
}
\value{A list with 7 elements:
\item{formula}{The model formula used in the function call.}
\item{model.dimension}{The model dimension under Poisson sampling.}
\item{status}{Description of how the algorithm terminated. For debugging purposes.}
\item{iterations}{Number of iterations the algorithm used.}.
\item{face}{The contingency table with an appended column denoting whether each cell is part of the facial set (denoted 0) or not (denoted 1). Those cells that are not part of the facial set are cells with an estimated cell mean of 0. These cells cause the non-existence of the MLE and should be removed from the table before the model is fit by the glm function. See the vignette for more details.}
\item{face.dimension}{The dimension of the facial set. Equivalent to the number of log-linear parameters that can be estimated.} 
\item{maxloglik}{The maximum of the log likelihood function up to a constant.} 
}
\references{Feinberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in log-linear models. Annals of Statistics. 40: 996-1023

Friedlander, M. (2016). Fitting log-linear models in sparse contingency tables using the eMLEloglin R package. Preprint. arXiv:1611.07505}
\author{Matthew Friedlander} 
\examples{
data(rochdale)
f <- facial_set (data = rochdale, 
       formula = freq ~ a*d + a*e + b*e + c*e + e*f + a*c*g + d*g + f*g + b*d*h)
N <- sum(rochdale$freq)
BIC <- f$maxloglik - f$model.dimension * log(N) / 2
cBIC <- f$maxloglik - f$face.dimension * log(N) / 2 
}
\keyword{htest}
\keyword{models}
