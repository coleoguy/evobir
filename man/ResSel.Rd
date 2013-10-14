\name{ResSel}
\alias{ResSel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Selection on Residuals}
\description{
This function takes measurements of multiple traits and performs a linear regression and identifies those records with the largest and smallest residual.  Originally it was written to perform a regression of horn size on body size allowing for high and low selection lines.}
\usage{
ResSel(data, traits, percent = 10, identifier = 1, model = "linear")
}
\arguments{
  \item{data}{this is a dataframe with subject identifiers and phenotypic trait values}
  \item{traits}{a numeric vector indicating the column containing the predictor and response variables in that order}
  \item{percent}{the percentage of highest and lowest residuals that should be identified}
  \item{identifier}{the column which contains the record numbers to identify individuals}
  \item{model}{currently this is not used}
}
\value{
This function returns a list
\item{high line}{ the ID numbers for the individuals selected for the high line}
\item{low line}{ the ID numbers of the individuals selected for the low line}
}
\references{\url{http://www.uta.edu/karyodb/evobiR/}}
\author{Heath Blackmon}
\examples{
data <- read.csv(file = system.file("horn.beetle.csv", package = "evobiR"))
ResSel(data = data, traits = c(2,3), percent = 15, identifier = 1, model = "linear")
}
\keyword{ hard selection }
\keyword{ residual }
