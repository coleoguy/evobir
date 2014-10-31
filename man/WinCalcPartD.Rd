\name{WinCalcPartD}
\alias{WinCalcPartD}
\title{Calculate Partitioned D-statistic in sliding windows}
\description{This functions calculate Partitioned D-statistic in windows.}
\usage{
WinCalcPartD(alignment = "alignment.fasta", win.size = 1000, , step.size=500,
boot = F, replicate = 1000)
}
\arguments{
  \item{alignment}{This is an alignment in fasta format}
  \item{win.size}{This is the size of the window used}
  \item{step.size}{This is the size of steps in the sliding window}
  \item{boot}{This indicates whether or not bootstrapping should be performed to estimate variance}
  \item{replicate}{Number of replicates to be used in estimating variance}
}
\details{ This function is just an extension of CalcPartD and calculates the partitioned D statistic for in windows. }
\value{
Returns a table with the number of each type of site, Z scores and p-values for each window in the genome
}
\author{Heath Blackmon}
\references{http://www.uta.edu/karyodb/evobiR/}

\examples{
WinCalcPartD(alignment = system.file("2.fasta", package = "evobiR"), 
          win.size=1000, step.size=500, boot = F, replicate=100)  }
\keyword{ introgression }
\keyword{ SNP }
\keyword{ ABBA }
\keyword{ sliding window }
\keyword{ D-statistic }
\keyword{ hybridization }
