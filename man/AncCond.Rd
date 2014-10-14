\name{AncCond}
\alias{AncCond}
\title{Calculate the mean of a continuous character at the origin of derived state of a binary charachter}
\description{This function uses stochastic mapping and ancestral state reconstruction to determine if the derived state of a binary trait originates when a continuous trait has an extreme value.}
\usage{
AncCond(trees, data, derived.state, iterations=1000)
}
\arguments{
  \item{trees}{tree(s) of class phylo or multiPhylo}
  \item{data}{a dataframe with rows matching the tips on the tree, and 2 columns the first should have the continuous trait value and the second the state of the binary character}
  \item{derived.state}{the derived condition for the binary trait}
  \item{iterations}{the number of iterations to be used in estimating significance}
}
\details{This function uses stochastic mapping and ancestral state reconstruction as implemented in phytools to determine if the derived state of a binary trait originates when a continuous trait has an extreme value.  This test assumes that the derived state of the binary character may lead to correlated selection in the continuous trait.  Because of this the ancestral state reconstruction of the continuous trait is based only on data from species that remain in the ancestral condition for the binary trait}
\value{
Returns a plot of the null distribution and the observed data as well as empirical p-value for the observed data.
}
\author{Heath Blackmon and Richard H. Adams}
\references{http://www.uta.edu/karyodb/evobiR/
\examples{
AncCond(trees, data, derived.state, iterations=100)  
}
\keyword{ comparative method }
\keyword{ discrete trait }
\keyword{ continuous trait }