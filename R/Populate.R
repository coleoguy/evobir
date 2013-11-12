Populate <- function(census, genome) {    # this will create an initial population of males and females 
  population <- list()
  for(i in 1:census){
    foo <- genome
    if(i <= {census / 2}){  # makes the first half males by deleting one X chromosome
      foo[c("x.2.s","x.2.d"), ] <- 0
      foo["xy.recom2", 1] <- 0
    }
    if(i > {census / 2}){         # makes the second half females by deleting the Y chromosome
      foo[c("y.1.d","y.1.s"), ] <- 0 
    }
    population[[i]]<-foo
  }
  return(population)
}