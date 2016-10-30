\name{getData}
\alias{getData}
\title{
Facilitates extraction of continuous data to match tips on a tree}
\description{
Supplied with a tree in phylo format and a csv file this function will use a column to find matching data for the tree.  In cases where more than one row in the csv file matches the tree tip it has the option of selecting a value from a random mathcing row, taking the mean of all matching rows, or the mode of matching rows.}
\usage{
getData(tree=tree, dat=dat, id.col=1, tr.col=c(3,5), method="R")
}
\arguments{
  \item{tree}{a phylogeny in phylo format}
  \item{dat}{the number of parameters in the model}
  \item{id.col}{the sample size.}
  \item{tr.col}{the sample size.}
  \item{method}{This is a string and must specify the method for handling multiple matches.  Options are "R" for a random match, "Me" for the mean of all matches, "Mo" for the mode of all matches.}
}
\details{
This function returns a list with the first object being a pruned version of the tree and the second being either a named numeric vector if only a single trait is given.  If more than one trait is selected then it will return a matrix with row names mathcing the tree and each column representing data for one of the target traits.}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}

\examples{
}
\keyword{ tree manipulation }
\keyword{ data munging }
