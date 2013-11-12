CrossoverMaps <- function(loci, genome.length){
  map.set <- list()
  for(i in 1:200){
    autosome.crossovers <- sort(sample(1:loci, sample(c(1, 2, 3), 1, prob = c(.25, .5, .25))))
    mp <- vector()
    if(length(autosome.crossovers) == 1){
      mp[1:autosome.crossovers[1]] <- 0
      mp[autosome.crossovers[1]:loci] <- 2
    }
    if(length(autosome.crossovers) == 2){
      mp[1:autosome.crossovers[1]] <- 0
      mp[autosome.crossovers[1]:autosome.crossovers[2]] <- 2
      mp[autosome.crossovers[2]:loci] <- 0
    }
    if(length(autosome.crossovers) == 3){
      mp[1:autosome.crossovers[1]] <- 0
      mp[autosome.crossovers[1]:autosome.crossovers[2]] <- 2
      mp[autosome.crossovers[2]:autosome.crossovers[3]] <- 0
      mp[autosome.crossovers[3]:loci] <- 2
    }
    if(loci < genome.length) mp[(loci+1):genome.length] <- 0
    map.set[[i]] <- mp
  }
  return(map.set)
}
    