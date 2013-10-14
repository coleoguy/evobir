\name{HetLevels}
\alias{HetLevels}
\title{
Calculates heterozygosity based on ambigous IUPAC nucleotide codes}
\description{
Some assemblers like AllPaths will use IUPAC codes to indicate polymorphic sites in an assembly.  This function uses these to estimate the level of heterozygosity along a sequence.}
\usage{
HetLevels(fasta)
}
\arguments{
  \item{fasta}{A fasta formatted file that may contain 1 or more sequences.}
}
\details{
Ambiguity codes of "r", "y", "s", "w", "k", "m","b", "d", "h", "v" are counted as polymorphic sites with an "n" are not counted since they are often indicative of scaffolding not polymorphism.}
\value{
A matrix is returned listing the length of each sequence the number of polymorphic sites and the number of "n" sites.
}
\references{
\url{http://www.uta.edu/karyodb/evobiR/}}
\author{
Heath Blackmon
}

\examples{
HetLevels(fasta = system.file("assembly.fasta", package = "evobiR"))
}
\keyword{ polymorphism }
\keyword{ snp }
