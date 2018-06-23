\name{HighLevelTree}
\alias{HighLevelTree}
\title{
  Prune a phylogeny to a higher taxonomic levels.
}
\description{
  This function takes a species level tree and prunes it so only a single representative of each clade, e.g. genus or family, is present.}
\usage{
  HighLevelTree(taxa.table, tree, cur.tips, new.tips)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxa.table}{a character matrix with columns containing the current tips of the tree and additional columns containing their assigned genera, family, tribe, etc.}
  \item{tree}{the current species level tree}
  \item{cur.tips}{a number indicating which column contains the current tip names}
  \item{new.tips}{a number indicating which column contains the names for the taxonomic level that you would like to prune the tree to.}
}
\value{
  returns and object of the class "phylo"
}
\references{
  \url{http://www.uta.edu/karyodb/evobiR/}}
\author{Heath Blackmon}
\seealso{ \code{\link{GetTaxonomy}}
}
\examples{
  data(species.tree)
  data(taxa.table)
  
  family.tree <- HighLevelTree(taxa.table, species.tree, 1, 2)
  genus.tree <- HighLevelTree(taxa.table, species.tree, 1, 3)
}
\keyword{ phylogeny }
\keyword{ comparative methods }