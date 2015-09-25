\name{ViewEvo}
\alias{ViewEvo}
\title{Learning Resources}
\description{
This uses the shiny app to produce interactive pages.
}
\usage{
ViewEvo(simulation)
}
\arguments{
  \item{simulation}{Text string indicating the application to run.  Currently options are "wf.model", "bm.model", "bd.model", "dist.model"}
}
\details{
The wf.model was implemented to illustrate to students the effects of genetic drift.  In particular the high likelihood of losing a beneficial allele when population size is finite.  

The bm.model was implemented to illustrate Brownian motion to students. It provides plots of multiple iterations of a trait evolving by brownian motion over the specified number of generations.  It also allows the student to plot the distribution of iterations at any selected point in time.  

The bd.model will plot 4 phylogenetic trees based on a birth death model with a single set of parameters.  This application was developed to illustrate the high variability of a birth death process as a generating model for phylogenies and the inherint difficulty in detecting differential diversification rates.  

Finally the dist.model was developed to help illustrate the relationship between common statistical distributions often used as priors and the way that parameters effect the density distribution.}
\value{This function returns an interactive webpage.}
\references{
\url{http://coleoguy.github.io/}

Wright-Fisher Simulator: \url{https://evobir.shinyapps.io/wf_model/}

Brownian Motion Simulator: \url{https://evobir.shinyapps.io/bm_model/}

Birth-death Simulator: \url{https://evobir.shinyapps.io/bd_model}

Statistical Distribution: \url{https://evobir.shinyapps.io/dist_model}
}
\author{Heath Blackmon}
\examples{
\dontrun{
ViewEvo("wf.model")
ViewEvo("bm.model")
ViewEvo("bd.model")
ViewEvo("dist.model")
}
}
\keyword{ pedagogical }
\keyword{ interactive }
\keyword{ teaching }
