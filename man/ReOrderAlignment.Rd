\name{ReOrderAlignment}
\alias{ReOrderAlignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Re-order sequences based on starting position}
\description{
This function re-orders aligned sequences with sequences based on the site of the alignment that they first have data. It also allows user to select reference sequences that will stay at the top of the alignmment regardless of starting position.  
}
\usage{
ReOrderAlignment(file, newfile, ref = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{the name of the fasta file which has sequences need to be re-ordered}
  \item{newfile}{name of the output file}
  \item{ref}{Either an integer or a characrter vector. Elements of the vector should match the sequences that should be kept at the top of the alignment.}
}
\value{A fasta file}
\references{http://coleoguy.github.io/}
\author{Sean Chien and Heath Blackmon}