\name{AncCond}
\alias{AncCond}
\title{Calculate the mean of a continuous character at transitions in a binary charachter}
\description{This function uses ancestral state estimations for a discrete character based on stochastic mapping under an Mk2 model and ancestral state estimates for a continuous trait under a Brownian motion model to determine if transitions in the binary trait coincide with extreme values of a continuous trait.}
\usage{
AncCond(tree, 
        data, 
        mc = 1000, 
        drop.state=NULL, 
        mat=c(0,2,1,0), 
        pi="equal", 
        n.tails = 1, 
        message = TRUE)
}
\arguments{
  \item{tree}{tree of class phylo}
  \item{data}{a dataframe with 3 columns. The first should match the taxa names in the tree, the second should have the continuous trait values and the third the states for the binary character. The binary trait should be coded as 1 and 2 if one is ancestral then it must be coded as 1.}
    \item{mc}{the number of iterations to use in building the null distribution.}
    \item{drop.state}{\code{NULL} or a numeric value of 1 or 2. If 1 or 2 are given then continuous data from species in the specified state will be dropped from the reconstruction of the continuous character.}
    \item{mat}{a vector describing the possible transition for the discrete trait. The default is equivelant to APE's \code{"ARD"} model \code{c(0,2,1,0)}, other options are \code{c(0,0,1,0)} which would allow only for transition from state 1 to state 2 or \code{c(0,1,1,0)} which would allow for transition in either direction but at equal rates.}
    \item{pi}{The probabilities for the binary trait at the root of the tree. The values possible are \code{"equal"}, \code{"estimated"}, or a numeric vector of length 2 with probabilities for each state}
    \item{n.tails}{numeric value of 1 or 2 to indicate whether a 1 or 2 tailed p-value should be calculated}
    \item{message}{Logical value if \code{TRUE} then status messages will be printed to the console}

}
\details{This function uses ancestral state estimates to determine if the transitions in the binary trait are associated with extreme values of the continuous trait.  This test can incorporate the possibility that the derived state of a binary character may lead to correlated selection in the continuous trait.  If this is desired then the drop.state argument should be used to specify the derived state of the binary character that should not be used in the ancestral state estimation of the continuous trait.}

\value{Returns a list of length 4:
  \item{OriginatingVals}{the mean value(s) for the continuous trait at the transition points of the binary character}
  \item{NTrans}{the number of transitions in the binary character}
  \item{NullDist}{the null distribution(s) produced by simulation}
  \item{pval}{pvalue}


}
\author{Nathan Anderson, Jeffery P. Demuth, Richard H. Adams, and Heath Blackmon}
\references{http://coleoguy.github.io/}
\examples{
\dontrun{
data(mite.trait)
data(trees.mite)
AncCond(trees[[1]], mite.trait) 
}
}
\keyword{ comparative method }
\keyword{ discrete trait }
\keyword{ continuous trait }