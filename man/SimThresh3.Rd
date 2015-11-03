\name{SimThresh3}
\alias{SimThresh3}
\title{
Simulates data for a 3 state discrete character evolving under a threshold model
\description{
This assumes that a 2D liability value determines which of 3 discrete states is observed.  The liabilities are allowed to evolve on an X and Y axis and then projected onto a cartesian plane.  This plane is divided into three regions and angles are defined measuring clockwise from with the origin as the vertex.  

330-90 = state 1
90-210 = state 2
210-330 = state 3

this approach allows for any set of probabilites for transitions between states to be represented.
}
\usage{
SimThresh3(tree, liabilities=F)}
\arguments{
  \item{tree}{a phylo object}
  \item{liabilities}{T or F to indicate whether underlying liabilities should be returned.}
}
\details{
If liabilities is set to TRUE then this function will return a list of length 3 containing the observed states the x values and y values of each tip in elements 1 through 3 respectively.
}
\references{
\url{http://coleoguy.github.io/}}
\author{
Heath Blackmon
}
\keyword{ comparative methods }
\keyword{ discrete trait evolution }
