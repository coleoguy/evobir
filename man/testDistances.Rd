\name{getDist}
\alias{getDist}
\title{test distances between changes}
\description{
This function takes a phylogeny and two binary traits (values should be coded as 1 and 2).  It then uses stochastic maps and an anova to determine whether any pair of trait transitions are closer together than expect if both traits are evolving neutrally.  A positive result suggests that there is a form of contingency between the traits.}
\usage{
testDistances(tree, trait1, trait2, n=100, model="ARD", verbose=F)
}
\arguments{
  \item{tree}{a phylo object}
  \item{trait1}{a named vector where the names match the tree and values are either 1 or 2.}
  \item{trait2}{a named vector where the names match the tree and values are either 1 or 2.}
  \item{n}{a named vector where the names match the tree and values are either 1 or 2.}
  \item{model}{describes the rate model for each of the binary traits.  Like \code{ACE} options are "ARD" and "SYM"}
  \item{verbose}{If TRUE then the distance estimates that are analyzed by ANOVA are returned.  If FALSE only the p-value of the ANOVA is returned}

}
\value{
returns either a p-value or a p-value and a list of distances.}

\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}
\keyword{ trait evolution }