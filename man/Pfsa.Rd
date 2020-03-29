\name{Pfsa}
\alias{Pfaa}
\alias{Pfss}
\title{Calculate the proportion of different classes of fusions as a fraction of all fusions}
\description{Thes three functions, Pfsa, Pfaa, and Pfss use sex chromosome systems and diploid autosome number to calculate the proportion of all fusions that are expected to be sex chromosome autosome fusions (Pfsa), autosome autosome fusions (Pfaa), and sex chromosome sex chromosome fusions (Pfss). }
\usage{
Pfsa(Da = 26, scs = "XY", mud = 0.5)
Pfaa(Da = 26, scs = "XY", mud = 0.5)
Pfss(Da = 26, scs = "XXY", mud = 0.5)
}
\arguments{
# Each of the three functions below contains three arguments
# Da: diploid count of autosomes
# scs: scs system possible values are
#      XO, XXO, XY, XYY, XXY, ZO, ZZO, ZW, ZWW, ZZW
# mud: proportion of fusions that originate in females
#      should be 0.5 unless other data is avaialable


  \item{Da}{the diploid autosome count}
  \item{scs}{a text string describing the sex chromosome system for instance: "XO", "XXO", "XY", "ZO"", or "ZWW"}
    \item{mud}{the proportion of fusions that originate in females.}

}
\details{This approach assumes that all chromosomes have equal probability of being involved in a fusion and that X and Y (or Z and W) chromosomes are not able to fuse. It will provide accurate results for any sex chromosome system that has any number of one sex chromosomes and either zero or one of the other. For instance it will work for XO, XY, XXY, XYYYY, but not for XXYY. It is applicable to both male and female heterogameic systems.}

\value{Returns a numeric vector of length 1 which contains the proportion of fusions expected to be of the specified type.}
\author{Nathan Anderson and Heath Blackmon}
\references{http://coleoguy.github.io/}
\examples{Pfsa(Da = 26, scs = "XY", mud=0.4)}
\keyword{ sex chromosome }
\keyword{ fusion }
