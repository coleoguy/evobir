\name{ReorderData}
\alias{ReorderData}
\title{
Reorders trait data to match the order of tips in a tree}
\description{
This function takes a vector, matrix, or dataframe and reorders the data to match the order of tips in a phylo object.}
\usage{
ReorderData(tree, data, taxa.names="row names")
}
\arguments{
  \item{tree}{a phylo object}
  \item{data}{a vector, matrix, dataframe set of taxa names as present in the tree and data must match.  If data is a vector it should be a named vector.  If the data is a matrix or dataframe the taxa names may be row names or present in a column.}
  \item{taxa.names}{If taxa names are present in a column the column number should be supplied.  If taxa names are the row names the argument can be set to "row names" (default setting).  If the data is being supplied in a vector this argument is not used.}
}
\details{
Returns data in the same format as supplied but reordered to match the order of tip labels.}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}
\keyword{ comparative analysis }
\keyword{ data munging }
