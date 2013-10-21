\name{PPSDiscrete}
\alias{PPSDiscrete}
\title{Create Simulated Datasets via PPS
}
\description{
This function performs posterior predictive simulations of discrete traits.  The function is written to work with the output of bayesian programs that produce a collection of rate matrix parameter estimates based on either one or a collection of trees.}
\usage{
PPSDiscrete(trees, MCMC, states, N = 2)
}
\arguments{
  \item{trees}{an object of class "multiPhylo" or "phylo" containing the trees used in generting the rate estimates}
  \item{MCMC}{this will normally be a log file that is brought into R with read.csv the columns for a three state character should be: tree, qAA, qBA, qCA, qAB, qBB, qCB, qAC, qBC, qCC.  If your analysis involves only a single tree then the tree column should be excluded.}
  \item{states}{a vector of root probabilities}
  \item{N}{the number of PPS datasets desired}
}
\value{
A matrix is returned with the rownames being the species names from the tree and each column containing a result of a single PPS.
}
\references{\url{http://www.uta.edu/karyodb/evobiR/}}
\author{Heath Blackmon}
\examples{
data(trees)
data(mcmc2)
data(mcmc3)
# 1 tree 100 q-mats 3 states
PPSDiscrete(trees[[1]], MCMC=mcmc3[,2:10], states=c(.5,.2,.3), N=2)
# 10 trees 100 q-mats 3 states
PPSDiscrete(trees, MCMC=mcmc3, states=c(.5,.2,.3), N=10)
# 10 trees 100 q-mats 2 states
PPSDiscrete(trees, MCMC=mcmc2, states=c(.5,.5), N=10)
}
\keyword{ PPS }
