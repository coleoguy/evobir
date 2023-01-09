\name{ReOrderAlignment}
\alias{ReOrderAlignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Re-order sequences by physical map}
\description{
This function re-orders aligned sequences based on physical map.  
}
\usage{
ReOrderAlignment(file, newfile, ref = (""))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{the name of the fasta file which has sequences need to be re-ordered}
  \item{newfile}{name of the output file}
  \item{ref}{Either an integer or a characrter vector. Elements of the vector should match the sequences that are wanted to be on the top.}
}
\value{A fasta file}
\references{http://coleoguy.github.io/}
\author{Sean Chien and Heath Blackmon}