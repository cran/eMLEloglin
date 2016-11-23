\name{eMLEloglin}
\Rdversion{1.0}
\title{Fitting log-linear models in sparse contingency tables.}
\description{Log-linear modeling is a popular method for the analysis of contingency table data. When the table is sparse, the data can fall on the boundary of the convex support, and we say that "the MLE does not exist" in the sense that some parameters cannot be estimated. However, an extended MLE always exists, and a subset of the original parameters will be estimable. The 'eMLEloglin' package determines which sampling zeros contribute to the non-existence of the MLE. These problematic zero cells can be removed from the contingency table and the model can then be fit using the glm() function.}
\details{
\tabular{ll}{
Package: \tab eMLEloglin\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-11-21\cr
License: \tab GPL-2\cr
}
The function facial_set identifies sampling zeros that contribute to the non-existence of the MLE.}
\author{Author: Matthew Friedlander
Maintainer: Matthew Friedlander <mattyf5@hotmail.com>
}
\references{Feinberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in log-linear models. Annals of Statistics, 40, 996-1023.}
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

