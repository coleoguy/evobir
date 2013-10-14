\name{CoalSim}
\alias{CoalSim}
\title{Coalescence Simulations}
\description{This function simulates the coalescence of linages through time.  At the begining of the simulation a number of unique lineages specified by the census size as we go back in time these coalesce until we have only a single lineage remaining.}
\usage{
CoalSim(census, lw, ln.col)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{census}{ the size of population to use in the simulation.}
  \item{lw}{ line width factor passed to the plotting function.}
  \item{ln.col}{ specifies the color to be used in plotting the coalesence history}
}
\value{
  Returns a plot showing the result of the coalescence simulation.
}
\references{http://www.uta.edu/karyodb/evobiR/}
\author{Heath Blackmon}
\examples{
CoalSim(census = 15, lw = 2, ln.col = 'blue')
}
\keyword{coalescence}
\keyword{population genetics}
