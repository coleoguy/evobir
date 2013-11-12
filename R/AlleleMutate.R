AlleleMutate <- function(x, y, sex.ant, loci.a, loci.s, DFE) {     # Apply a point mutation to an individual    # x = population    y = result of function RadDose
  foobish <- x[[y[[2]]]]  #this pulls out the genome we are mutating
  bias <- sample(c(0, 1, 2), size = 1, prob = c(sex.ant[1], sex.ant[2], sex.ant[3])) # 0 = no sex bias; 1 = male bias; 2 = female bias
  locus.draw <- sample(1:{loci.a + loci.s}, size = 1)
  ## first mutations withou sex bias
  if(bias == 0){
    if(locus.draw <= loci.a){      #autosome mutation
      bar <- sample(c(1, 3), 1)
      foobish[c(bar, bar + 1), locus.draw] <- sample(DFE, 1)
    }else{
      locus.draw <- locus.draw - loci.a
      if(y[[1]] == "Male"){    # male genome sex chromosomes
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
      foobish[sample(c(1, 3), 1), locus.draw] <- sample(DFE, 1)
    }else{
      locus.draw <- locus.draw - loci.a
      if(y[[1]] == "Male"){    # male genome sex chromosomes
        foobish[sample(c(5, 9), 1), locus.draw] <- sample(DFE, 1)
      }else{
        foobish[sample(c(5, 7), 1), locus.draw] <- sample(DFE, 1)
      }
    }
  }
  if(bias == 2){       ##   now mutations with female expression
    if(locus.draw <= loci.a){      
      foobish[sample(c(2, 4), 1), locus.draw] <- sample(DFE, 1)
    }else{
      locus.draw <- locus.draw - loci.a
      if(y[[1]] == "Male"){    # female expressed sex chromosomes
        foobish[sample(c(6, 10), 1), locus.draw] <- sample(DFE, 1)
      }else{
        foobish[sample(c(6, 8), 1), locus.draw] <- sample(DFE, 1)
      }
    }
  }
  x[[y[[2]]]] <-foobish  #this puts the mutated genome back in the population
  return(x)
}
