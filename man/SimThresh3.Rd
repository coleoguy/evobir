\name{SimThresh3}
\alias{SimThresh3}
\title{
3 state threshold model
}
\description{This is an extension of Felsenstein threshold model for a 3 state unordered trait.  This assumes that a 2 dimensional liability value determines which of 3 discrete states is observed.  The liabilities are allowed to evolve on an X and Y axis and then projected onto a cartesian plane.  This plane is divided into three sectors that meet at the origin.  The discrete trait observed is determined by which sector the X,Y values fall in: 330 to 90 degrees is state 1, 90 to 210 degrees is state 2, and 210 to 330 is state 3.  This approach allows for any set of probabilites for transitions between states to be represented.
}
\usage{SimThresh3(tree, liabilities=F, sig2, root)}
\arguments{
  \item{tree}{a phylo object}
  \item{liabilities}{TRUE or FALSE to indicate whether underlying liabilities should be returned.}
  \item{sig2}{instantaneous variance for Brownian motion model}
  \item{root}{numeric vector of length two that specifies the starting liabilities}
}

\details{
If liabilities is set to TRUE then this function will return a list of length 3 containing the observed states the x values and y values of each tip in elements 1 through 3 respectively. If set to false then a named numeric vector containing tip states (1 to 3) is returned.
}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}
\keyword{ comparative methods }
\keyword{ discrete trait evolution }
