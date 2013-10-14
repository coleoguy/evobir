GenomeSym <- function(model              = "fragileY",    # fragileY, canonical, distance, full
generations        = 100,
census             = 50,            # total census size
gametes            = 4,             # number of gametes drawn from each individual
loci.a             = 100,           # loci on autsomes
loci.s             = 50,            # loci on sex chromosomes
prob.point.mut     = .0010,         # probability that an new allele will arrise
DFE                = (c(rnorm(1000,mean=.9,sd=.06),rgamma(1000, rate=2, shape=.4))),    # hist(DFE,main="DFE",xlab="Fitness",xlim=c(0,1.5),breaks=40)
sex.bias.mut       = .85,           # propotion of mutations to arrise in males
haplosuff.mut      = 1,             # mutations in haplosufficiency will be drawn from a uniform distribution between 0 and N
crossovers         = 2,             # on average N crossovers per chromosome per meiosis in autosomes
recom.mu           = 10,             # maximum number of loci N that can quit recombining in a single mutational step
anneuploi.mut      = c(0, .0002),   # probability of a Y aneuploidy mutation drawn from a uniform dist between X and Y
distance.mut       = c(0, .000002), # prob of evolving dist pairing uniform dist from X to Y
achiasmat.mut      = c(0, .00002),  # prob of evolving achiasmatic meiosis in males unirom from X to Y
sex.ant            = c(.5, .3, .2), # prob of unbiased, male biased, female biased
fragile.fact       = 20,            # this is the factor by which we increase aneuploidy based on par reduction
reporting          ="None"){         # PAR, XY.fitness, all.loci

##############################################
## ALL INTERNAL FUNCTIONS TILL AROUND LINE 360
##############################################

  
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

DesignLife <- function(loci.a) {    # this will create a genome file that our population will start with
    genome <- matrix(1, 14, loci.a)
    # naming convention for the rows is auto.1.s is autosome copy 1 fitness in sires
    # x1.d is the fitness of x chromosome copy 1 in a dam
    rownames(genome) <- c("auto.1.s", "auto.1.d", "auto.2.s", "auto.2.d", 
                          "x.1.s",    "x.1.d",    "x.2.s",    "x.2.d", 
                          "y.1.s", "y.1.d",                                                         # we have to have the y.1.d because recombination can turn thes into X alleles
                          "xy.recom",                                                               # thinking about how recombination ceases this should be the union of the parents
                          "haplosufficiency.x1","haplosufficiency.x2","haplosufficiency.y1")        # haplo-sufficiency should be inherited with an X or Y
    genome[12:14, ] <- 0                                                                            # we begin with no haplosufficiency
    return(genome)
  }
  
Populate <- function(census) {    # this will create an initial population of males and females 
    population <- list()
    for(i in 1:census){
      foo <- genome
      if(i <= (census / 2)){  # makes the first half males by deleting one X chromosome
        foo["x.2.s", ] <- 0
        foo["x.2.d", ] <- 0
        foo["haplosufficiency.x2",] <- 0
        foo["xy.recom", 1] <- 0 # makes the sex determining locus non recombining
        # this is really just a necesity for the way that
        # I allow recombination to evolve obviously having only
        # not recombining doesnt really do anything.
      }
      if(i > (census / 2)){         # makes the second half females by deleting the Y chromosome
        foo["y.1.s", ] <- 0 
        foo["y.1.d", ] <- 0 
        foo["haplosufficiency.y1",]<-0 
        foo["xy.recom", 1] <- 0     # see above 4 lines
      }
      population[[i]]<-foo
    }
    return(population)
  }
  
sex <- function(x) {   # this function checks the sex of a genome
    if(sum(x[7, ]) == 0){
      if(sum(x[8, ]) == 0){ # we dont check the Y here that way we can have XO males
        foo <- "Male"
      }
    }
    if(sum(x[7, ]) != 0){
      if(sum(x[8, ]) != 0){
        if(sum(x[9, ]) == 0){
          foo <- "Female"
        }
      }
    }
    return(foo)
  }

RadDose <- function(x) {   # this will pick who to assign a mutation to  # x is the matrix produced by the function Ranking
    x <- as.data.frame(x)
    for(i in 1:nrow(x)){
      if(x[i, 1] == "Male"){                         # these two if statements prob of mutations that should occur in males vs females
        x[i, 3] <- sex.bias.mut
      }
      if(x[i, 1] == "Female"){
        x[i, 3] <- 1 - sex.bias.mut
      }
    }
    foo <- sample(1:census, replace = F, prob = x[, 3])   # here we do the sampling
    bar <- list()
    bar[1] <- as.character(x[foo[1], 1])
    bar[2] <- foo[1]
    return(bar)                                    # here we just return a list with the sex and the index number
  }
  
CeaseRecomb <- function(x, y) {   # Mutations that extend the SLR(sex limited region)   # x = population    y = result of function RadDose
    foobish <- x[[y[[2]]]]                         # this pulls out the genome we are mutating
    inv.size <- sample(1:recom.mu, 1)
    inv.start <- which.max(foobish[11, ] == 1) - 1
    foobish[11, (inv.start):(inv.start + inv.size)] <- 0
    x[[y[[2]]]] <- foobish
    return(x)
  }
  
ChangeHaploSuff <- function(x, y) {   # Mutations that extend the SLR(sex limited region)   # x = population    y = result of function RadDose
    newhap.val <- sample(runif(20, min = 0, max = 1), 1)
    foobish <- x[[y[[2]]]]  #this pulls out the genome we are mutating
    which.locus <- sample(1:loci.s, 1)
    if(sex(foobish) == "Male"){
      which.hap <- sample(c(12, 14), 1)
      foobish[which.hap, which.locus] <- newhap.val
    }else{
      which.hap <- sample(c(12, 13), 1)
      foobish[which.hap, which.locus] <- newhap.val
    }
    foobish -> x[[y[[2]]]]
    return(x)
  }
  
AlleleMutate <- function(x, y) {     # Apply a point mutation to an individual    # x = population    y = result of function RadDose
    foobish <- x[[y[[2]]]]  #this pulls out the genome we are mutating
    bias <- sample(c(0, 1, 2), size = 1, prob = c(sex.ant[1], sex.ant[2], sex.ant[3])) # 0 = no sex bias; 1 = male bias; 2 = female bias
    locus.draw <- sample(1:(loci.a + loci.s), size = 1)
    ## first mutations withou sex bias
    if(bias == 0){
      if(locus.draw <= loci.a){      #autosome mutation
        bar <- sample(c(1, 3), 1)
        foobish[c(bar, bar + 1), locus.draw] <- sample(DFE, 1)
      }else{
        locus.draw <- locus.draw - loci.a
        if(sport[[1]] == "Male"){    # male genome sex chromosomes
          bar <- sample(c(5, 9), 1)
          foobish[c(bar, bar + 1),locus.draw] <- sample(DFE, 1)
        }else{
          bar <- sample(c(5, 7), 1)    # female genome sex chromosomes
          foobish[c(bar, bar + 1), locus.draw] <- sample(DFE, 1)
        }
      }
    }  
    if(bias == 1){      ##   now mutations with male sex bias
      if(locus.draw <= loci.a){      
        bar <- sample(c(1, 3), 1)
        foobish[bar, locus.draw] <- sample(DFE, 1)
      }else{
        locus.draw <- locus.draw - loci.a
        if(y[[1]] == "Male"){    # male genome sex chromosomes
          bar <- sample(c(5, 9), 1)
          foobish[bar, locus.draw] <- sample(DFE, 1)
        }else{
          bar <- sample(c(5, 7), 1)    # female genome sex chromosomes
          foobish[bar, locus.draw] <- sample(DFE, 1)
        }
      }
    }
    if(bias == 2){       ##   now mutations with female expression
      if(locus.draw <= loci.a){      
        bar <- sample(c(2, 4), 1)
        foobish[bar, locus.draw] <- sample(DFE, 1)
      }else{
        locus.draw <- locus.draw - loci.a
        if(y[[1]] == "Male"){    # female expressed sex chromosomes
          bar <- sample(c(6, 10), 1)
          foobish[bar, locus.draw] <- sample(DFE, 1)
        }else{
          bar <- sample(c(6, 8), 1)    # female expressed sex chromosomes
          foobish[bar, locus.draw] <- sample(DFE, 1)
        }
      }
    }
    foobish -> x[[y[[2]]]]  #this puts the mutated genome back in the population
    return(x)
  }
  
Ranking <- function(x) {  # this will evaluate the fitness of each individual    # x is the list of genomes in the population
  fitness <- function(x){      # this function calculates fitness
    if(sex(x) == "Male"){
      a.fitt <- prod(colMeans(x[c(1, 3), ]))       # here is the autosome fitness     
      s1 <- colMeans(x[c(5, 9), 1:loci.s ])       # means for sex chromosomes
      s2 <- apply(x[c(5, 9), 1:loci.s], 2, max)   # max for sex chromosomes
      hap <- x[c(12, 14), 1:loci.s]               # this sets up a two 2xloci.s matrix
      hap.ind <- apply(x[c(5, 9), 1:loci.s], 2, which.max) # this holds the index we need
      hap.val <- vector()
      for(k in 1:loci.s){                         
        hap.val[k] <- hap[hap.ind[k],k]
      }
      s.fitt <- prod(s1 + (s2 - s1) * hap.val)    # this is the SC fitness with haplosufficiency
      tot.fit <- s.fitt * a.fitt
    }
    if(sex(x) == "Female"){
      a.fitt <- prod(colMeans(x[c(2, 4), ]))       # here is the autosome fitness     
      s1 <- colMeans(x[c(6, 8), 1:loci.s ])       # means for sex chromosomes
      s2 <- apply(x[c(6, 8), 1:loci.s], 2, max)   # max for sex chromosomes
      hap <- x[c(12, 13), 1:loci.s]               # this sets up a two 2xloci.s matrix
      hap.ind <- apply(x[c(6, 8), 1:loci.s], 2, which.max) # this holds the index we need
      hap.val <- vector()
      for(k in 1:loci.s){                         
        hap.val[k] <- hap[hap.ind[k],k]
      }
      s.fitt <- prod(s1 + (s2 - s1) * hap.val)    # this is the SC fitness with haplosufficiency
      tot.fit <- s.fitt * a.fitt
    }
    return(tot.fit)
  }
  
  hotness <- matrix(, census, 2)
  colnames(hotness) <-  c("sex", "fitness")
  hotness[,1] <- unlist(lapply(x, sex))
  foo <- unlist(lapply(x, fitness))
  # here we need to divide foo values for average for that sex
  hotness[,2] <- as.numeric(foo)
  males <- hotness[, 1] == 'Male'
  male.foo <- mean(as.numeric(hotness[males, 2]))
  hotness[males,2] <- as.numeric(hotness[males,2])/male.foo
  female.foo <- mean(as.numeric(hotness[!males, 2]))
  hotness[!males,2] <- as.numeric(hotness[!males,2])/female.foo
  return(hotness)
}

GametoGenesis <- function(x) {        # Makes eggs and sperm       # x = population
    #internal functions
    Crossover <- function(x, dist, loci){                      # given current chrom (0, 2) , dist between rows of matrix and length of chrom calls crossovers with mean equal to spec.
      if(x == 0){                                               # I need this if so it can have a high probability of staying in current state...lame
        foo <- sample(c(0, dist), 1, c(loci / crossovers, 1), replace = T)  # this is the opportunity to cross over
      }else{
        foo <- sample(c(dist, 0), 1, c(loci/crossovers, 1), replace = T)  # this is the opportunity to cross over
      }
      return(foo)
    } 
    
    
    
    
    Spermatogenesis <- function(x, n){         #x is the fathers genome     n is the number of sperm to produce should be even
      ## FIRST WE TAKE CARE OF THE AUTOSOMES
      gams <- list()
      for(j in 1:(n/2)){                                
        mp <- 0                                             # this is the switch variable it can be 0 or 2 and a switch causes crossover
        gam1 <- gam2 <- matrix(0, 8, loci.a)
        for(i in 1:loci.a){                                 # this will run along our chromosomes and allow recombination to occur
          gam1[c(1, 2), i] <- x[c(1 + mp, 2 + mp), i]       # this builds the recombinant autosomal chromosome one locus at a time
          gam2[c(1, 2), i] <- x[c(3 - mp, 4 - mp), i]               # this builds the other resulting recombinant
          mp <- Crossover(mp, 2, loci.a)
        }
        ## NOW LETS DO THE SEX CHROMOSOMES gam1 = X bearing      gam2 = Y bearing
        mp <- 0                                            # this is the switch variable we start out at 0 no crossover  
        if(model != "distance"){
          for(i in 1:loci.a){                              # this will run along our chromosomes and allow recombination to occur        
            gam1[c(3, 4), i] <- x[c(5 + mp, 6 + mp), i]    # this builds the recombinant sex chromosome one locus at a time
            gam2[c(5,6),i] <- x[c(9-mp,10-mp),i]            # this builds the other resulting recombinant
            gam1[7,] <- gam2[7,] <- x[11,]  
            if(mp == 0){
              gam1[8,] <- x[12,]
              gam2[8,] <- x[14,]
            }else{
              gam1[8,] <- x[14,]
              gam2[8,] <- x[12,]
            }      
            if(x[11,i] == 1){                     ############# this if checks to see whether we have shut down recombination in only calls possible crossover when we have the ability to recombine
              mp <- Crossover(mp, 4, loci.s)
            }
          }
        }else{
          gam1[c(3, 4, 7, 8), ] <- x[c(5, 6, 11, 12), ]
          gam2[c(5, 6, 7, 8), ] <- x[c(9, 10, 11, 14), ]
        }
        gams[[j]] <- gam1
        gams[[(n / 2) + j]] <- gam2
      }
      return(gams)
    }
    
    
    
    
    Oogenesis <- function(x, n){    #x is the mothers genome     n is the number of sperm to produce should be even
      ## FIRST WE TAKE CARE OF THE AUTOSOMES
      gams <- list()
      for(j in 1:(n / 2)){                                
        mp <- 0                                             # this is the switch variable it can be 0 or 2 and a switch causes crossover
        gam1 <- gam2 <- matrix(, 6, loci.a)
        for(i in 1:loci.a){                                 # this will run along our chromosomes and allow recombination to occur
          gam1[c(1, 2), i] <- x[c(1 + mp, 2 + mp), i]               # this builds the recombinant autosomal chromosome one locus at a time
          gam2[c(1, 2), i] <- x[c(3 - mp, 4 - mp), i]               # this builds the other resulting recombinant
          mp <- Crossover(mp, 2, loci.a)
        }
        ## NOW LETS DO THE SEX CHROMOSOMES gam1 = X bearing      gam2 = Y bearing
        mp <- 0                                            # this is the switch variable we start out at 0 no crossover  
        for(i in 1:loci.a){                              # this will run along our chromosomes and allow recombination to occur        
          gam1[c(3, 4), i] <- x[c(5 + mp, 6 + mp), i]            # this builds the recombinant autosomal chromosome one locus at a time
          gam2[c(3, 4), i] <- x[c(7 - mp, 8 - mp), i]            # this builds the other resulting recombinant
          gam1[5, ] <- gam2[5, ] <- x[11, ]  
          if(mp == 0){
            gam1[6, ] <- x[12, ]
            gam2[6, ] <- x[13, ]
          }else{
            gam1[6, ] <- x[13, ]
            gam2[6, ] <- x[12, ]
          }      
          mp <- Crossover(mp, 2, loci.s)
        }
        gams[[j]] <- gam1
        gams[[(n / 2) + j]] <- gam2
      }
      return(gams)
    } 
    
    
    
    GametePool <- list()
    for(i in 1:census){
      if(evaluation[i, 1] == "Male"){
        foo <- Spermatogenesis(x[[i]], gametes)
      }else{
        foo <- Oogenesis(x[[i]], gametes)
      }
      GametePool[[i]] <- foo
    }
    bar<-list()
    k<- 1
    for(i in 1:census){
      for(j in 1:gametes){
        bar[[k]] <- GametePool[[i]][[j]]
        k <- k + 1
      }
    }
    return(bar)
  }
  
NextGen <- function(x, y) {  # Assembles the next generation with selection   # x = gamete.pool    y = evaluation
    new.pop <- list()
    # SELECTION ACTS NOW BY EFFECTING THE PROBABILITY OF LEAVING OFFSPRING
    foo <- which(y[, 1] == "Male")
    foo2 <- which(y[, 1] == "Female")
    dads<-moms<-vector()
    for(i in 1:census){
      dads[i] <- (sample(foo, 1, prob = y[foo, 2], replace = T) * gametes) - sample(0:(gametes - 1), 1, replace = T)   # select male gametes based on fathers fitness
      moms[i] <- (sample(foo2, 1, prob = y[foo2, 2], replace = T) * gametes) - sample(0:(gametes - 1), 1, replace =T ) # select female gametes based on mohters fitness
    }
    # NOW LETS ASSEMBLE OUR NEW GENOMES
    genome <- matrix(0,14,loci.a)
    rownames(genome) <- c("auto.1.s", "auto.1.d", "auto.2.s", "auto.2.d", "x.1.s",    "x.1.d",    "x.2.s",    "x.2.d", 
                          "y.1.s", "y.1.d", "xy.recom", "haplosufficiency.x1","haplosufficiency.x2","haplosufficiency.y1")
    for(i in 1:census){
      baz <- genome
      # lets inherit XY recombination as the more restricted of the gametes
      baz[11, max(c(which.max(x[[moms[i]]][5, ] == 1), which.max(x[[dads[i]]][7, ] == 1))):loci.s] <- 1
      # and the rest of fertilization
      baz[c(3, 4, 5, 6, 12), ] <- x[[moms[i]]][c(1, 2, 3, 4, 6), ]
      baz[c(1, 2, 7, 8, 9, 10, 14), ]  <- x[[dads[i]]][c(1, 2, 3, 4, 5, 6, 8), ]
      new.pop[[i]]<-baz
    }
    return(new.pop)
  }

#################################################                   
## HERE BEGINS THE ACTUAL CODE FOR RUNNING THE SIM.
#################################################
                   
  genome <- DesignLife(loci.a)            # This creates the matrix for an individual genome
  population <- Populate(census)          # This creates a population 1:1 sex ration of size=census
  if(reporting != "None") {               # activate any in simulation visualization
                     if(reporting == 'all.loci'){
                       y.chrom <- matrix(, loci.s, generations)
                       layout(matrix(c(1,2),2,1))
                       plot(0.1,0.1, xlim = c(0,generations), ylim = c(.95,1.05), cex = .05, main = 'Locus Mean Fitness', xlab = 'Gen.', ylab = 'Abs. Fitness')
                       p1 <- par(no.readonly=T)
                       plot(0.1,0.1, xlim = c(0,generations), ylim = c(.5,2), cex = .05, main = 'Locus Mean Fitness', xlab = 'Gen.', ylab = 'Abs. Fitness')
                       p2 <- par(no.readonly=T)
                       loci.col <- rainbow(100)
                     }
                     if(reporting == 'xy'){
                       y.chrom <- matrix(, loci.s, generations)
                       layout(matrix(c(1,2),1,2))
                       plot(0.1,0.1, xlim = c(0,generations), ylim = c(.5,2), cex = .05, main = 'X Chromosome', xlab = 'Gen.', ylab = 'Abs. Fitness')
                       text(0, .61, 'FEMALE', col='red', pos=4, font=2)
                       text(0, .50, 'MALE', col='blue', pos=4, font=2)
                       p1 <- par(no.readonly=T)
                       plot(0.1,0.1, xlim = c(0,generations), ylim = c(.5,2), cex = .05, main = 'Y Chromosome', xlab = 'Gen.', ylab = 'Abs. Fitness')
                       p2 <- par(no.readonly=T)
                       Sys.sleep(.1)
                       
                     }
                     if(reporting == "PAR")
                       plot(0,0,ylim=c(0,50),xlim=c(0,generations),col="white",pch=19,main="XY Recombination", ylab="Number of Loci Isolated", xlab="Generations")
                     if(reporting == "XY.fitness")
                       plot(0,0,ylim=c(0,50),xlim=c(0,generations),col="green",pch=19,main="XY Recombination", ylab="Number of Loci Isolated", xlab="Generations")
                   }
    for(k in 1:generations){                ## THIS BEGINS THE GENERATION SIMULATION PROCESS
      cat("Generation: ",k," has been born.\n",sep="")
      evaluation <- Ranking(population)       #This sexes and evaluates the fitness of all individuals
      if(reporting != "None") {        # update any in simulation visualization
        if(reporting == "PAR") {
                         foo <- vector()
                         for(j in 1:census){
                           foo[j] <- which.min(population[[j]][11,] == 0) - 1
                         }
                         points(k,mean(foo),col="darkgreen",pch=19,cex=.5)
                         Sys.sleep(.1)
                       }
        if(reporting == 'all.loci'){
                         par(p1)
                         par(mfg=c(1,1))
                         foo <- matrix(,length(population),loci.a)
                         for(i in 1:length(population)){
                           foo[i,] <- population[[i]][1,]
                         }
                         points(rep(k,loci.a),colMeans(foo), cex = .2, col = loci.col)
                         par(p2)
                         par(mfg=c(2,1))
                         points(rep(k,loci.a),colMeans(foo), cex = .2, col = loci.col)
                         Sys.sleep(.25)
                       }
        if(reporting == 'xy'){
                         sires <- which(evaluation[,1] == 'Male')
                         dams <- which(evaluation[,1] == 'Female')
                         par(p1)
                         par(mfg=c(1,1))      
                         # lets do average X in a female first
                         fooxd <- fooxs <- vector()
                         for(i in 1:length(dams)){
                           fooxs[i] <- mean(population[[dams[i]]][5, 1:loci.s])
                           fooxd[i] <- mean(population[[dams[i]]][6, 1:loci.s])
                         }
                         points(rep(k, length(fooxd)), fooxd, cex = .2, col = 'red')      
                         points(rep(k, length(fooxs)), fooxs, cex = .2, col = 'blue')      
                         fooy <- vector()
                         for(i in 1:length(sires)){
                           fooy[i] <- mean(population[[sires[i]]][9, 1:loci.s])
                         }
                         par(p2)
                         par(mfg=c(1,2))
                         points(rep(k, length(fooy)), fooy, cex = .2, col = 'blue')
                       }
        #Sys.sleep(.3)
      }
      # Premiotic Mutations -----    
      if(sample(1:100,1) > 1) {
        sport <- RadDose(evaluation)                       # This picks an individual from the population to be mutated
        population <- AlleleMutate(population, sport)      # This performs a mutation ceating a new allele
      }
      if(sample(1:100,1) > 50) {
        sport <- RadDose(evaluation)                       # This picks an individual from the population to be mutated
        population <- CeaseRecomb(population, sport)       # Mutations that extend the SLR(sex limited region)
      }
      if(sample(1:100,1) > 95) {
        sport <- RadDose(evaluation)                       # This picks an individual from the population to be mutated
        population <- ChangeHaploSuff(population, sport)   # Mutations that change halposufficiency
      }
      # Meiosis and Fertilization -----    
      gamete.pool <- GametoGenesis(population)           # This creates the number of gametes specified
      population <- NextGen(gamete.pool, evaluation)     # This creates the population for generation i+1
    }
return(population)
}

                   