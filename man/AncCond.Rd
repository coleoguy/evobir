\name{AncCond}
\alias{AncCond}
\title{Calculate the mean of a continuous character at the origin of derived state of a binary charachter}
\description{This function uses stochastic mapping and ancestral state reconstruction to determine if the derived state of a binary trait originates when a continuous trait has an extreme value.}
\usage{
AncCond(trees, data)
}
\arguments{
  \item{trees}{tree of class phylo}
  \item{data}{a dataframe with 3 columns. The first should match the taxa names in the tree, the second should have the continuous trait values and the third the states for the binary character.  The binary character should be coded as 0 and 1. With 0 as the ancestral state and 1 as the derived state.}
    \item{iter}{the number of iterations to use in builing the null distribution. Defaults to 1000.}

}
\details{This function uses ancestral state reconstruction to determine if the derived state of a binary trait originates when a continuous trait has an extreme value.  This test assumes that the derived state of the binary character may lead to correlated selection in the continuous trait.  Because of this the ancestral state reconstruction of the continuous trait is based only on data from species that remain in the ancestral condition for the binary trait}
\value{Returns a list of length 4:
  \item{OriginatingVals}{the mean value for the continuous trait at the origin of the derived state of the binary character}
  \item{NTrans}{the number of origins for the derived state of the binary character}
  \item{NullDist}{the null distribution produced by simulation}
  \item{pval}{pvalue}


}
\author{Heath Blackmon and Richard H. Adams}
\references{http://coleoguy.github.io/}
\examples{
\dontrun{
data(mite.trait)
data(trees.mite)
AncCond(trees, mite.trait) 
}
}
\keyword{ comparative method }
\keyword{ discrete trait }
\keyword{ continuous trait }