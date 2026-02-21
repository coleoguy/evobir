\name{SuperMatrix}
\alias{SuperMatrix}
\title{
Create a supermatrix from multiple gene alignments}
\description{
Combines all alignments in a folder into a single supermatrix. Missing data is inserted for taxa that are absent from individual alignments. Partition boundaries are tracked so the user knows which columns correspond to which input alignment.}
\usage{
SuperMatrix(missing = "-", prefix = "concatenated", save = T,
input = "", format.in = "f", format.out = "f", concatenate = T)
}
\arguments{
  \item{missing}{the character to use when no data is available for a taxon. Defaults to \code{"-"}.}
  \item{prefix}{prefix for the resulting supermatrix output files. Defaults to \code{"concatenated"}.}
  \item{save}{logical; if \code{TRUE} then supermatrix and partition files will be saved to the working directory. Defaults to \code{TRUE}.}
  \item{input}{a regular expression determining which files will be read in. For example, \code{".fasta"} will read in all files containing ".fasta" in their names. Default is blank and will result in all files in the working directory being read in.}
  \item{format.in}{a character string specifying the format of the alignments to be read in. Passed to \code{read.dna} in \pkg{ape}: \code{"interleaved"}, \code{"sequential"}, \code{"clustal"}, or \code{"fasta"}, or abbreviations of these.}
  \item{format.out}{a character string specifying the format for the supermatrix to be saved to. Passed to \code{write.dna} in \pkg{ape}: \code{"interleaved"}, \code{"sequential"}, \code{"clustal"}, or \code{"fasta"}, or abbreviations of these.}
  \item{concatenate}{logical; when \code{TRUE} sequences are concatenated into a single alignment. When \code{FALSE} sequences are saved as individual files but are expanded to include all taxa in the combined dataset. Defaults to \code{TRUE}.}
}
\details{
This function reads all alignments in the working directory (or those matching the \code{input} pattern) and constructs a single supermatrix that includes all taxa present in any of the files and inserts missing symbols for taxa that are missing sequences for some loci.}
\value{
A list with two elements:
  \item{partitions}{a dataframe with columns for partition name, start position, and stop position recording the alignment boundaries of each input file in the combined supermatrix}
  \item{supermatrix}{a character matrix of the combined supermatrix}

If \code{save = TRUE}, both files are also saved to the working directory.
}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}

\examples{
\dontrun{
SuperMatrix(missing = "N", prefix = "DATASET2", save = TRUE,
            format.in = "f", format.out = "f")
}
}
\keyword{ concatenation }
\keyword{ alignment }
\keyword{ supermatrix }
