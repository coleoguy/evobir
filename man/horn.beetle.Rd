\name{horn.beetle.csv}
\alias{horn.beetle.csv}
\docType{data}
\title{Gnatocerus cornutus measurements}
\description{A csv file containing measurements of horn and body size for the broad-horned flour beetle \emph{Gnatocerus cornutus}. The file has three columns: individual ID, body size, and horn size. Used in the \code{\link{ResSel}} function examples for demonstrating residual selection analysis.}
\examples{
data <- read.csv(file = system.file("horn.beetle.csv", package = "evobiR"))
head(data)
}
\keyword{datasets}
