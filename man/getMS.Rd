\name{getMS}
\alias{getMS}
\title{
creates an MS string}
\description{
Supplied with a phylo object it returns a character string describing a tree for ms.}
\usage{
getEJ(tree, samps, report, n.site)
}
\arguments{
  \item{tree}{a phylo format tree}
  \item{samps}{number of samples desired}
  \item{report}{either "T", or "s" for output of gene trees or segregating sites.}
  \item{n.site}{number of segregating sites desired.}
}
\details{
Returns a character vector}
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
