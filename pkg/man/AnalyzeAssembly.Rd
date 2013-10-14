\name{AnalyzeAssembly}
\alias{AnalyzeAssembly}
\title{Analyze Genome Assembly}
\description{This function reports basic statistics for a genome assembly.}
\usage{
AnalyzeAssembly(genome, max_N = 25, plot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genome}{ a list of vectors with each element being a single string of the class "SeqFastadna".}
  \item{max_N}{ Maximum number of consecutive N symbols.  Scaffolds will be broken into contigs 
  when this number is exceeded.}
  \item{plot }{ When True an accumulation plot will be returned as well as the statistics}
}

\details{
If a standard FASTA file is read in with the function read.fasta from the package seqinr the argument as.string should set to TRUE.  The genome should also be all lower case which is the default setting for read.fasta.}
\value{
  A dataframe with the following rows:
  
  Number of Scaffolds\cr
  Assembly Size Based on Scaffolds\cr
  Number of Scaffolds over 1MB\cr
  N50 Scaffold Size\cr
  Number of Contigs\cr
  Assembly Size Based on Contigs\cr
  N50 Contig Size\cr
  Minimum Contig Size\cr
  Percent GC\cr
}
\references{http://www.uta.edu/karyodb/evobiR/}
\author{Heath Blackmon}
\examples{
## just a small simulated genome
data(genome)
## calculate summary statistics for the genome
AnalyzeAssembly(genome = genome, max_N = 25, plot = TRUE)
}
\keyword{assembly}
\keyword{genome}
