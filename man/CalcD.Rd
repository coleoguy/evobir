\name{CalcD}
\alias{CalcD}
\alias{CalcPopD}
\title{Calculate Patterson's D-statistic}
\description{These functions calculate Patterson's D-statistic to compare the frequencies of discordant SNP genealogies.  These tests assume equal substitution rates and unlinked loci, D-statistics significantly different from 0 suggest that introgression has occurred.}
\usage{
CalcD(alignment = "alignment.fasta", sig.test = "N", block.size = 1000, replicate = 1000)

CalcPopD(alignment = "alignment.fasta")

}
\arguments{
  \item{alignment}{This is an alignment in fasta format.  Sequences should be in the order: P1, P2, P3, Outgroup.}
  \item{sig.test}{This indicates whether or if to test for significance. Options are "B" bootstrap, "J" jackknife, or "N" none.}
  \item{block.size}{The number of sites to be dropped in the jackknife approach}
  \item{replicate}{Number of replicates to be used in estimating variance}
}
\details{ The functions CalcD and CalcPopD are implementations of the algorithm described in Durand et al. 2011.  Significance of the D-stat can be calculated either through bootstrapping or jackknifing.  Bootstrapping is appropriate for datasets where SNPs are unlinked for instance unmapped RADSeq data.  Jackknifing is the appropriate approach when SNPs are potentially in linkage for instance gene alignments or genome alignments.}
\value{
Returns the number of each type of site, Z scores and p-values
}
\author{Heath Blackmon}
\references{http://coleoguy.github.io/

Durand, Eric Y., et al. Testing for ancient admixture between closely related populations. Molecular biology and evolution 28.8 (2011): 2239-2252.

Eaton, D. A. R., and R. H. Ree. 2013. Inferring phylogeny and introgression using RADseq data: An example from flowering plants (Pedicularis: Orobanchaceae). Syst. Biol. 62:689-706}

\examples{
CalcD(alignment = system.file("1.fasta", package = "evobiR"), sig.test = "N")

CalcPopD(alignment = system.file("3.fasta", package = "evobiR"))
  }
\keyword{ introgression }
\keyword{ SNP }
\keyword{ ABBA }
\keyword{ D-statistic }
\keyword{ hybridization }
