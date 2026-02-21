\name{CalcD}
\alias{CalcD}
\alias{CalcPopD}
\title{Calculate Patterson's D-statistic}
\description{These functions calculate Patterson's D-statistic to compare the frequencies of discordant SNP genealogies. These tests assume equal substitution rates and unlinked loci; D-statistics significantly different from 0 suggest that introgression has occurred. \code{CalcD} uses single sequences per taxon while \code{CalcPopD} uses population-level data with multiple sequences per population.}
\usage{
CalcD(alignment = "alignment.fasta", sig.test = "N", ambig="D",
block.size = 1000, replicate = 1000, align.format='fasta')

CalcPopD(alignment = "alignment.fasta", sig.test = "N", ambig="D",
block.size = 1000, replicate = 1000, align.format='fasta')

}
\arguments{
  \item{alignment}{path to an alignment file. Sequences should be in the order: P1, P2, P3, Outgroup. In the case of \code{CalcPopD}, sequences from each population should have identical names (see file 3.fasta for an example).}
  \item{sig.test}{indicates whether and how to test for significance. Options are \code{"B"} bootstrap, \code{"J"} jackknife, or \code{"N"} none.}
  \item{ambig}{flag indicating how to deal with sequence ambiguities. Options are \code{"D"} drop all ambiguous loci, \code{"R"} resolve each biallelic ambiguity randomly, or \code{"I"} ignore ambiguity and perform analysis without checking sequences. If the argument \code{"R"} is chosen there is a degree of stochasticity in the analysis and the user should run the analysis more than once. It is also wise to compare this result to setting the argument to \code{"D"}.}
  \item{block.size}{the number of sites to be dropped in the jackknife approach}
  \item{replicate}{number of replicates to be used in estimating variance}
  \item{align.format}{a character string specifying the format of the alignment file: \code{"mase"}, \code{"clustal"}, \code{"phylip"}, \code{"fasta"}, or \code{"msf"}}
}
\details{The functions \code{CalcD} and \code{CalcPopD} are implementations of the algorithm described in Durand et al. 2011. Significance of the D-statistic can be calculated either through bootstrapping or jackknifing. Bootstrapping is appropriate for datasets where SNPs are unlinked, for instance unmapped RADSeq data. Jackknifing is the appropriate approach when SNPs are potentially in linkage, for instance gene alignments or genome alignments.}
\value{
\code{CalcD} returns the D-statistic as a numeric value.

\code{CalcPopD} returns a list with the following elements:
  \item{d.stat}{the D-statistic}
  \item{pval}{the p-value (only present when \code{sig.test} is not \code{"N"})}
  \item{align.length}{the number of sites in the alignment}
  \item{useful.sites}{the number of sites with ABBA or BABA patterns}
  \item{seg.sites}{the number of ABBA/BABA sites still segregating in at least one population}

Both functions also print summary statistics to the console.
}
\author{Heath Blackmon}
\references{http://coleoguy.github.io/

Durand, Eric Y., et al. Testing for ancient admixture between closely related populations. Molecular biology and evolution 28.8 (2011): 2239-2252.

Eaton, D. A. R., and R. H. Ree. 2013. Inferring phylogeny and introgression using RADseq data: An example from flowering plants (Pedicularis: Orobanchaceae). Syst. Biol. 62:689-706}

\examples{
CalcD(alignment = system.file("1.fasta", package = "evobiR"), sig.test = "N")

CalcPopD(alignment = system.file("3.fasta", package = "evobiR"), sig.test = "N")
  }
\keyword{ introgression }
\keyword{ SNP }
\keyword{ ABBA }
\keyword{ D-statistic }
\keyword{ hybridization }
