Fitness <- function(x, loci.s){      # this function calculates fitness x = population
  if(sex(x) == "Male"){
    a.fitt <- prod(colMeans(x[c(1, 3), ]))      # here is the autosome fitness     
    s1 <- colMeans(x[c(5, 9), 1:loci.s ])       # means for sex chromosomes
    s2 <- apply(x[c(5, 9), 1:loci.s], 2, max)   # max for sex chromosomes
    hap <- x[c(13, 14), 1:loci.s]               # this sets up a two 2xloci.s matrix
    hap.ind <- apply(x[c(5, 9), 1:loci.s], 2, which.max) # this holds the index we need
    hap.val <- vector()
    for(k in 1:loci.s){                         
      hap.val[k] <- hap[hap.ind[k],k]
    }
    s.fitt <- prod(s1 + {s2 - s1} * hap.val)    # this is the SC fitness with haplosufficiency
    tot.fit <- s.fitt * a.fitt
  }
  if(sex(x) == "Female"){
    a.fitt <- prod(colMeans(x[c(2, 4), ]))      # here is the autosome fitness     
    s1 <- colMeans(x[c(6, 8), 1:loci.s ])       # means for sex chromosomes
    s2 <- apply(x[c(6, 8), 1:loci.s], 2, max)   # max for sex chromosomes
    hap <- x[c(13, 14), 1:loci.s]               # this sets up a two 2xloci.s matrix
    hap.ind <- apply(x[c(6, 8), 1:loci.s], 2, which.max) # this holds the index we need
    hap.val <- vector()
    for(k in 1:loci.s){                         
      hap.val[k] <- hap[hap.ind[k],k]
    }
    s.fitt <- prod(s1 + {s2 - s1} * hap.val)    # this is the SC fitness with haplosufficiency
    tot.fit <- s.fitt * a.fitt
  }
  return(tot.fit)
}