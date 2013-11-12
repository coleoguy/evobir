Spermatogenesis <- function(x, n, loci.a, model, loci.s, crossovers, auto.mp.set){         #x is the fathers genome     n is the number of sperm to produce should be even
  ## FIRST WE TAKE CARE OF THE AUTOSOMES
  ## gam1 female gam2 male
  mp <- auto.mp.set[[sample(1:length(auto.mp.set), 1)]]
  gams <- list()
  gam1 <- gam2 <- matrix(0, 8, loci.a)                               
  for(i in 1:loci.a){                                               # this will run along our chromosomes and allow recombination to occur
    gam1[c(1, 2), i] <- x[c(1 + mp[i], 2 + mp[i]), i]               # this builds the recombinant autosomal chromosome one locus at a time
    gam2[c(1, 2), i] <- x[c(3 - mp[i], 4 - mp[i]), i]               # this builds the other resulting recombinant
  }
  sex.chr.strt <- which.max(x[12,] == 1)            # this finds the first region that can have a crossover
  if(sex.chr.strt > 48){
      gam1[c(3, 4, 7, 8), ] <- x[c(5, 6, 11, 13), ]                 # this builds the recombinant autosomal chromosome one locus at a time
      gam2[c(5, 6, 7, 8), ] <- x[c(9, 10, 12, 14), ]                # this builds the other resulting recombinant
    }else{
      mp.s <- mp/2                                                  # gam1 will be female gam2 will be male
      for(i in 1:{sex.chr.strt - 1}){                               # this will run along our chromosomes and allow recombination to occur        
        gam1[c(3, 4, 7, 8), i] <- x[c(5, 6, 11, 13), i]             # this builds the recombinant autosomal chromosome one locus at a time
        gam2[c(5, 6, 7, 8), i] <- x[c(9, 10, 12, 14), i]            # this builds the other resulting recombinant
      }
      for(i in sex.chr.strt:loci.s){                                # this will run along our chromosomes and allow recombination to occur        
        gam1[c(3, 4, 7, 8), i] <- x[c(5 + {2*mp[i]}, 6 + {2*mp[i]}, 11 + mp.s[i], 13 + mp.s[i]), i]            # this builds the recombinant autosomal chromosome one locus at a time
        gam2[c(5, 6, 7, 8), i] <- x[c(9 - {2*mp[i]}, 10 - {2*mp[i]}, 12 - mp.s[i], 14 - mp.s[i]), i]            # this builds the other resulting recombinant
      }
    }
  gams[[1]] <- gam1
  gams[[2]] <- gam2
  return(gams)
}