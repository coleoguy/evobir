Ranking <- function(x, census, loci.s) {  # this will evaluate the fitness of each individual    # x is the list of genomes in the population
  hotness <- matrix(, census, 2)
  colnames(hotness) <-  c("sex", "fitness")
  hotness[,1] <- unlist(lapply(x, sex))
  foo <- unlist(lapply(x, Fitness, loci.s=loci.s))
  # here we need to divide foo values for average for that sex
  hotness[,2] <- as.numeric(foo)
  males <- hotness[, 1] == 'Male'
  male.foo <- mean(as.numeric(hotness[males, 2]))
  hotness[males,2] <- as.numeric(hotness[males,2])/male.foo
  female.foo <- mean(as.numeric(hotness[!males, 2]))
  hotness[!males,2] <- as.numeric(hotness[!males,2])/female.foo
  return(hotness)
}