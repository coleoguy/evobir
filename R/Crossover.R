Crossover <- function(x, dist, loci, crossovers){                      # given current chrom (0, 2) , dist between rows of matrix and length of chrom calls crossovers with mean equal to spec.
  if(x == 0){                                               # I need this if so it can have a high probability of staying in current state...lame
    foo <- sample(c(0, dist), 1, c(loci / crossovers, 1), replace = T)  # this is the opportunity to cross over
  }else{
    foo <- sample(c(dist, 0), 1, c(loci/crossovers, 1), replace = T)  # this is the opportunity to cross over
  }
  return(foo)
} 
