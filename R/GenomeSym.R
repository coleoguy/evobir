source("AlleleMutate.R")
source('CeaseRecomb.R')
source('ChangeHaploSuff.R')
source('DesignLife.R')
source('GametoGenesis.R')
source('Mode.R')
source('NextGen.R')
source('Oogenesis.R')
source('Populate.R')
source('RadDose.R')
source('Ranking.R')
source('sex.R')
source('Spermatogenesis.R')
source('Crossover.R')
source('Fitness.R')
source('CrossoverMaps.R')
source('Monitor.R')

GenomeSymC <- function(model              = "canonical",    # fragileY, canonical, distance, full
generations        = 100,
census             = 50,            # total census size
gametes            = 2,             # number of gametes drawn from each individual
loci.a             = 50,            # loci on autsomes
loci.s             = 25,            # loci on sex chromosomes
point.mut          = 10,            # number of point mutations per generation
DFE                = c(rnorm(1000,mean=.9,sd=.06),rgamma(1000, rate=2, shape=.4)),    # hist(DFE,main="DFE",xlab="Fitness",xlim=c(0,1.5),breaks=40)
sex.bias.mut       = .85,           # propotion of mutations to arrise in males
haplosuff.mut      = 1,             # mutations in haplosufficiency will be drawn from a uniform distribution between 0 and N
haplo.rate         = 1,             # number of haplosufficiency mutations per generation                       
crossovers         = 2,             # on average N crossovers per chromosome per meiosis in autosomes
                                    # currently crossover is restricted to being 2 to change this need to update code in oogenesis.R and spermatogenesis.R
recom.mu           = 4,             # maximum number of loci N that can quit recombining in a single mutational step
rate.recom.mut     = 10,            # how many recombination mutations per generation
anneuploi.mut      = c(0, .0002),   # probability of a Y aneuploidy mutation drawn from a uniform dist between X and Y
distance.mut       = c(0, .000002), # prob of evolving dist pairing uniform dist from X to Y
achiasmat.mut      = c(0, .00002),  # prob of evolving achiasmatic meiosis in males unirom from X to Y
sex.ant            = c(.5, .3, .2), # prob of unbiased, male biased, female biased
fragile.fact       = 20,            # this is the factor by which we increase aneuploidy based on par reduction
report.freq        = 10,            # how frequently to monitor the population                       
report.style       = "population",   # the data to collect on the population: population, y-chrom, PAR-active, PAR-save                       
seed.val           = 1){            # to help tracking errors
  
  set.seed(seed.val)
  counter <- 0
  
  DFE[DFE>1.1] <- .95  # this is a stop gap I seem to have way to many really high fitness mutations
  Snowden <- list()                   
  genome       <- DesignLife(loci.a)                                              # This creates the matrix for an individual genome
  population   <- Populate(census, genome)                                        # This creates a population 1:1 sex ration of size=census
  auto.mp.set  <- CrossoverMaps(loci.a, loci.a)
  sex.mp.set   <- CrossoverMaps(loci.s, loci.a)
  evaluation <- Ranking(population, census, loci.s)                             # This sexes and evaluates the fitness of all individuals
  
  for(k in 1:generations){                                                        ## THIS BEGINS THE GENERATION SIMULATION PROCESS

    cat(".")
    if(k/25 == round(k/25)) cat(k,"\n")

    # Premiotic Mutations -----    
    for(i in 1:point.mut){
      population <- AlleleMutate(population, RadDose(evaluation, sex.bias.mut), 
                                 sex.ant, loci.a, loci.s, DFE)                    # This performs a mutation ceating a new allele
    }
    for(i in 1:rate.recom.mut){
      population <- CeaseRecomb(population, recom.mu, loci.s, evaluation)         # Mutations that extend the SLR(sex limited region)
    }
    for(i in 1:haplo.rate){
      population <- ChangeHaploSuff(population, 
                                    RadDose(evaluation, sex.bias.mut), 
                                    haplosuff.mut, loci.s)                        # Mutations that change halposufficiency
    }

    #Meiosis and Fertilization -----    
    gamete.pool <- GametoGenesis(population, crossovers, 
                                 census, evaluation, 
                                 gametes, loci.a, model, 
                                 loci.s, auto.mp.set, sex.mp.set)                 # This creates the number of gametes specified
    population <- NextGen(gamete.pool, evaluation, census, gametes, loci.a)       # This creates the population for generation i+1

    evaluation <- Ranking(population, census, loci.s)                             # This sexes and evaluates the fitness of all individuals
    
    if(round(k/report.freq) == k/report.freq){
      counter <- counter+1
      if(counter == 1) rep.cols <- rainbow(round(generations/report.freq))
      Snowden[[counter]] <- Monitor(population, report.style, report.freq, k, loci.s, generations, evaluation, counter, rep.cols)
    }
  }
return(Snowden)
}