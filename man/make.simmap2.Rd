\name{make.simmap2}
\alias{make.simmap2}
\title{Modified stochastic mapping which is resistant to model failure}
\description{This function is an extension of the \code{make.simmap} function of \pkg{phytools} which allows users to monitor rejections at specific edges during the simulation proces and set an upper limit on the number of rejections permitted per edge.}
\usage{
make.simmap2(tree, x, model="SYM", nsim=1 ,
             monitor = FALSE, rejmax = NULL, rejint = 1000000, 
             ...)
}
\arguments{
  \item{tree}{see \code{\link{make.simmap}}}
  \item{x}{see \code{\link{make.simmap}}}
  \item{model}{see \code{\link{make.simmap}}}
  \item{nsim}{see \code{\link{make.simmap}}}
    \item{monitor}{\code{boolean} describing whether or not to print number of rejections per edge to console. Defaults to \code{FALSE}}
  \item{rejmax}{\code{int} greater than one giving the maximum number of rejections permitted before an edge is skipped. Defaults to \code{NULL} (no upper limit on rejections)}
  \item{rejint}{\code{int} giving the interval after which rejection number is to be printed to console}
  \item{...}{optional arguments. see \code{\link{make.simmap}}}
}
\value{A object of class \code{simmap} or \code{multiSimmap}, see \code{\link{make.simmap}}. 

In addition to the states present in the model, an additional state \code{fail} in the \code{maps} and \code{mapped.edge} elements is assigned to edges which are skipped due to exceeding the rejection limit.
}
\seealso{
  \code{\link{make.simmap}}
}
\author{Matthew Marano, Maximos Chin, and Heath Blackmon}
\references{http://coleoguy.github.io/}
\keyword{comparative phylogenetics, discrete trait, stochastic mapping}
