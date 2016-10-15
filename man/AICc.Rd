\name{AICc}
\alias{AICc}
\title{
Computes an AICc score}
\description{
Supplied with a log likelihood, and the number of model parameters, and sample size calculates the small sample size version of the AIC score.}
\usage{
AICc(loglik, K, N)
}
\arguments{
  \item{loglik}{log likelihood.}
  \item{K}{the number of parameters in the model}
  \item{N}{the sample size.}
}
\details{
Returns an AICc score.}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}

\examples{
AICc(-32, 3, 100)
}
\keyword{ model comparison }
\keyword{ AICc }
