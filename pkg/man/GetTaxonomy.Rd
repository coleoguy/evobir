\name{GetTaxonomy}
\alias{GetTaxonomy}
\title{Create a taxonomy table}
\description{This function is an extension of the tax_name function from the package taxize.  The input is a tree and the function will return a table with a description of the taxonomy of all species on the tree}
\usage{GetTaxonomy(tree, database = "ncbi")}
\arguments{
  \item{tree}{an object of class "phylo"}
  \item{database}{character string "ncbi or "itis" indicating the database that should be used}
}
\value{
  this function saves a csv file named: tree.taxonomy.csv
}
\references{http://www.uta.edu/karyodb/evobiR/}
\author{Heath Blackmon}
\examples{
# data(hym.tree)
# GetTaxonomy(tree = hym.tree, database = 'ncbi')
}
\keyword{taxonomy}