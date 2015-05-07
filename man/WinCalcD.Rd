\name{WinCalcD}
\alias{WinCalcD}
\title{Calculate Patterson's D-statistic in sliding windows}
\description{This functions calculate Patterson's D-statistic in windows.}
\usage{
WinCalcD(alignment = "alignment.fasta", win.size = 500, , step.size=50,
boot = F, replicate = 1000)
}
\arguments{
  \item{alignment}{This is an alignment in fasta format}
  \item{win.size}{This is the size of the window used}
  \item{step.size}{This is the size of steps in the sliding window}
  \item{boot}{This indicates whether or not bootstrapping should be performed to estimate variance}
  \item{replicate}{Number of replicates to be used in estimating variance}
}
\details{ This function is just an extension of CalcD and calculates D statistic for windows. }
\value{
Returns a table with the number of each type of site, Z scores and p-values for each window in the genome
}
\author{Heath Blackmon}
\references{http://www.uta.edu/karyodb/evobiR/

Durand, Eric Y., et al. Testing for ancient admixture between closely related populations. Molecular biology and evolution 28.8 (2011): 2239-2252.

Eaton, D. A. R., and R. H. Ree. 2013. Inferring phylogeny and introgression using RADseq data: An example from flowering plants (Pedicularis: Orobanchaceae). Syst. Biol. 62:689-706}

\examples{
WinCalcD(alignment = system.file("1.fasta", package = "evobiR"), win.size=100, step.size=50, boot = TRUE, replicate=10)
  }
\keyword{ introgression }
\keyword{ SNP }
\keyword{ ABBA }
\keyword{ sliding window }
\keyword{ D-statistic }
\keyword{ hybridization }
