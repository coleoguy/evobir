Oogenesis <- function(x, n, loci.a, model, loci.s, crossovers, auto.mp.set, sex.mp.set){    #x is the mothers genome     n is the number of eggs to produce should be even
  # first we get a vector that represents the crossover events
  mp <- auto.mp.set[[sample(1:length(auto.mp.set), 1)]]
  gams <- list()
  gam1 <- gam2 <- matrix(1, 6, loci.a)
  for(i in 1:loci.a){                                               # this will run along our chromosomes and allow recombination to occur
    gam1[c(1, 2), i] <- x[c(1 + mp[i], 2 + mp[i]), i]               # this builds the recombinant autosomal chromosome one locus at a time
    gam2[c(1, 2), i] <- x[c(3 - mp[i], 4 - mp[i]), i]               # this builds the other resulting recombinant
  }
  #####THAT SHOULD TAKE CARE OF THE AUTOSOMES NOW WE NEED TO DO THE SEX CHROMOSOEMS
  mp <- sex.mp.set[[sample(1:length(sex.mp.set), 1)]]
  mp.s <- mp/2    
    for(i in 1:loci.a){                              # this will run along our chromosomes and allow recombination to occur        
      gam1[c(3, 4, 6), i] <- x[c(5 + mp[i], 6 + mp[i], 13 + mp.s[i]), i]            # this builds the recombinant autosomal chromosome one locus at a time
      gam2[c(3, 4, 6), i] <- x[c(7 - mp[i], 8 - mp[i], 14 - mp.s[i]), i]            # this builds the other resulting recombinant
     }
  # row five is left as straight ones since X chromosomes always recombine in females
  gams[[1]] <- gam1
  gams[[2]] <- gam2
  return(gams)
}
