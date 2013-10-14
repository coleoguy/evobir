\name{DriftSim}
\alias{DriftSim}
\alias{SelSim}
\title{Simulate Allele Frequencies}
\description{These functions simulate allele frequencies through time}
\usage{
SelSim(census, initial.freq, selection.coef, iter, generations)

DriftSim(census, initial.freq, iter, generations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{census}{ The population size to be used for the simulation.}
  \item{initial.freq}{ The proportion of individuals that should have the selected mutation at the begining of the simulation.}
  \item{selection.coef}{ The selection coefficient.  1 is a neutral mutation while values greater or less than 1 will be beneficial or deleterious respectively.}
  \item{iter}{ The number of times that the simulation should be repeated}
  \item{generations}{ How many generations each simulation should last.}
}
\value{
  These functions were created to help in teaching population genetics.  Each of these functions plot the frequency of a new mutation at each generation.  Each iteration of the simulation is plotted in a new color and the functions also reports the number of iterations that resulted in fixation, loss or intermediate frequencies for the mutation.  The function DriftSim is just a special case of SelSim where the selection coefficient is equal to 1.
}
\references{http://www.uta.edu/karyodb/evobiR/}
\author{Heath Blackmon}
\examples{
set.seed(5)
SelSim(census = 200, initial.freq = .5, selection.coef = 1.01, iter = 10, generations = 200)

set.seed(5)
DriftSim(census = 200, initial.freq = .5, iter = 25, generations = 200)
}
\keyword{population genetics}
\keyword{genetic drift}
\keyword{selection}