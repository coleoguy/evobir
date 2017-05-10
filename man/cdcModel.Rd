\name{cdcModel}
\alias{cdcModel}
\title{Fits a model of change dependent character evolution.}
\description{This function implements a hidden state in one binary trait and allows users to test whether transitions in a second binary trait occur primarily right after transitions in the first.}
\usage{
cdcModel(x, y, tree, model, initial.vals=NULL)
}
\arguments{
  \item{x}{first binary trait, this is the trait with a hidden state}
  \item{y}{second binary trait that may be correlated with the first}
  \item{tree}{tree of class phylo}
  \item{model}{a numeric vector either 333, 343, or 344 indicating the types of rates in character 2 (y).}
  \item{initial.vals}{the initial values to be used for find.mle. if not given they will all be set to .1}
}
\details{I'll explain this more someday}
\value{A list of class fit.mle

  - par estimats

  - lnLik the log likelihood of the ML point

  - counts the number of function evaluations

  - code 0 means the optimizere believes it found a solution

  - func the liklihood function used

  - method the optimisation method used
}
\author{Heath Blackmon}
\references{http://coleoguy.github.io/}
\examples{
\dontrun{
load(url("http://coleoguy.github.io/comparative/pagel.data.sets.RData"))
tree <- set.phy[[1]]
x <- unrepburst1[[1]]
y <- unrepburst2[[1]]
result <- cdcModel(x, y, tree, model=344)
}
}
\keyword{ comparative method }
\keyword{ discrete trait }
