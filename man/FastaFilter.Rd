\name{FastaFilter}
\alias{FastaFilter}
\title{
Filter fasta files}
\description{
Function to remove duplicate sequences.}
\usage{
FastaFilter(f.in, folder = F, prefix="pruned")
}
\arguments{
  \item{f.in}{Iif a single file is being processed then this is the filename including the relative path.}
  \item{folder}{if TRUE then all files in the current working directory will be processed if FALSE then a single file from f.in will be processed}
  \item{prefix}{A character string specifying the prefix to be added to the file names after filtering.}
}
\details{
Function removes duplicate sequences from the same taxon.  Will keep only the longest file for each taxon.  Can process a single file or all files in a folder.}
\value{
Saves a fasta formatted file that has been reduced to unique sequences.
}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}

\examples{

\keyword{ alignment DNA fasta }
