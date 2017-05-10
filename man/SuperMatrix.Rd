\name{SuperMatrix}
\alias{SuperMatrix}
\title{
creates a supermatrix from multiple gene alignments}
\description{
combines all alignments in a folder into a single supermatrix}
\usage{
SuperMatrix(missing = "-", prefix = "concatenated", save = T, 
format.in = "f", format.out = "f")
}
\arguments{
  \item{missing}{the character to use when no data is available for a taxa}
  \item{prefix}{prefix for the resulting supermatrix}
  \item{save}{if True then supermatrix and partition file will be saved}
  \item{format.in}{A character string specifying the format of the alignments to be read in. The argument is passed to read.dna in APE: "interleaved", "sequential", "clustal", or "fasta", or abbreviations of these are available.}
  \item{format.out}{A character string specifying the format for the supermatrix to be saved to. The argument is passed to write.dna in APE: "interleaved", "sequential", "clustal", or "fasta", or abbreviations of these are available.}  
  
  
}
\details{
This function reads all alignments in the working directory and constructs a single supermatrix that includes all taxa present in any of the files and inserts missing symbols for taxa that are missing sequences for some loci.}
\value{
A list with two elements is returned. The first element contains partition data that records the alignment positions of each input alignment file in the combined supermatrix.  The second element is a dataframe version of the supermatrix. If the argument save is set to TRUE then both of these files are also saved to the working directory.
}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}

\examples{
\dontrun{
SuperMatrix(missing = "N", prefix = "DATASET2", save = T, format = "f")
}
}
\keyword{ concatention alignment DNA fasta }
\keyword{ SuperMatrix, concatenation, alignment }
