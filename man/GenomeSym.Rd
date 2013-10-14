\name{GenomeSym}
\alias{GenomeSym}
\title{Genome Simulation}
\description{THIS FUNCTION IS STILL IN THE DEVELOPMENT PHASE AND RESULTS FROM IT SHOULD NOT BE TRUSTED WITHOUT CONTACTING ME FIRST!!!

This function is designed to simulate the evolution of a diploid organism with chromosomal sex determination.  It allows for the evolution of sexual antagonism, recombination supression, anneuploidy, haplosufficiency.  }
\usage{
GenomeSym(model = "fragileY", generations = 100, census = 50,  gametes = 4,
loci.a  = 100, loci.s = 50, prob.point.mut = .0010, DFE = (c(rnorm(1000, mean = 0.9, sd = 0.06), rgamma(1000, rate = 2, 
          shape = 0.4))), sex.bias.mut = .85, haplosuff.mut = 1, crossovers = 2, 
recom.mu = 10, anneuploi.mut = c(0, .0002), distance.mut = c(0, .000002), 
achiasmat.mut = c(0, .00002), sex.ant = c(.5, .3, .2), fragile.fact = 20, reporting = "None") }
\arguments{
  \item{model}{ the type of simulation to be run}
  \item{generations}{ integer indicating the number of generations that the simulation should be run}
  \item{census}{ integer indicating total census size}
  \item{gametes}{ integer indicating number of gametes drawn from each individual}
  \item{loci.a}{ integer indicating loci on autsomes}
  \item{loci.s}{ integer indicating loci on sex chromosomes}
  \item{prob.point.mut}{ probability that an new allele will arrise}
  \item{DFE}{ a vector of fitness effect for mutations}
  \item{sex.bias.mut}{ propotion of mutations to arrise in males}
  \item{haplosuff.mut}{ a value of less than 1 from which mutations in haplosufficiency will be drawn from a uniform distribution between 0 and the specified value.}
  \item{crossovers}{ the average number of crossovers per chromosome per meiosis in autosomes}
  \item{recom.mu}{ maximum number of loci N that can quit recombining in a single mutational step}
  \item{anneuploi.mut}{ probability of a Y aneuploidy mutation drawn from a uniform dist between X and Y}
  \item{distance.mut}{ prob of evolving dist pairing uniform dist from X to Y}
  \item{achiasmat.mut}{ prob of evolving achiasmatic meiosis in males unirom from X to Y}
  \item{sex.ant}{ numerical vector must sum to 1 indicating the probability of unbiased, male biased, female biased mutations}
  \item{fragile.fact}{ this is the factor by which we increase aneuploidy based on par reduction}
  \item{reporting}{ a character string indicating the desired type of insimulation reporting: 'PAR', 'XY.fitness', 'all.loci'.}
}
\details{
This function simulates a genome as a series of numerical vectors each chromosome in a genome has vectors to represent those genes fitness when in a female or when in a male additionally if they are sex chromosomes they have vectors describing the recombination characteristics and haplosufficiency data.  All of these genome characteristics are allowed to evolve during the simulation through a variety of mutational mechanisms.  The genomes perform in silico gametogenesis with recombination.  The probability of gametes being used to create the next generation is in proportion to their parents relative fitness.

default values for the function are as follows:
  
model              = "fragileY"

generations        = 100

census             = 50 

gametes            = 4   

loci.a             = 100    

loci.s             = 50     

prob.point.mut     = .3    

DFE                = (c(rnorm(1000,mean=.9,sd=.06),rgamma(1000, rate=2, shape=.4)))   

sex.bias.mut       = .85    

haplosuff.mut      = 1     

crossovers         = 2     

recom.mu           = 10      

anneuploi.mut      = c(0, .0002)   

distance.mut       = c(0, .000002) 

achiasmat.mut      = c(0, .00002)  

sex.ant            = c(.5, .3, .2) 

fragile.fact       = 20         

reporting          ="None"         
}
\references{http://www.uta.edu/karyodb/evobiR/}

\author{Heath Blackmon}
\examples{
## just a small simulated genome
Genome <- GenomeSym(census = 20, generations = 10, reporting = 'all.loci')
}
\keyword{assembly}
\keyword{genome}
