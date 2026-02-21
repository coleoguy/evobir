\name{ReOrderAlignment}
\alias{ReOrderAlignment}
\title{Re-order sequences based on starting position}
\description{
This function re-orders aligned sequences based on the site of the alignment where they first have data (i.e., the first non-gap position). It also allows the user to select reference sequences that will stay at the top of the alignment regardless of starting position. This is useful for visual inspection of alignments where sequences begin at different positions.
}
\usage{
ReOrderAlignment(file, newfile, ref = NULL)
}
\arguments{
  \item{file}{path to the fasta file containing sequences to be re-ordered}
  \item{newfile}{path for the output fasta file}
  \item{ref}{either an integer vector of sequence indices or a character vector of sequence names. Elements should identify the sequences that should be kept at the top of the alignment. Defaults to \code{NULL} (no reference sequences).}
}
\value{Writes a re-ordered fasta file to \code{newfile}. Called for its side effect.}
\references{http://coleoguy.github.io/}
\author{Sean Chien and Heath Blackmon}
\examples{
\donttest{
ReOrderAlignment(file = system.file("1.fasta", package = "evobiR"),
                 newfile = tempfile(fileext = ".fasta"))
}
}
