\name{SuperMatrix}
\alias{SuperMatrix}
\title{
creates a supermatrix from multiple gene alignments}
\description{
combines all alignments in a folder into a single supermatrix}
\usage{
superMatrix(missing = "-", format = "f",
            file = "conc.alignment.fasta")
}
\arguments{
  \item{missing}{the character to use when no data is available for a taxa}
  \item{format}{the format of the input }
  \item{N}{the sample size.}
}
\details{
Returns an AICc score.}
\references{
\url{http://www.uta.edu/karyodb/evobiR/}}
\author{
Heath Blackmon
}

\examples{
AICc(-32, 3, 100)
}
\keyword{ model comparison }
\keyword{ AICc }
