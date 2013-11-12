GametoGenesis <- function(x, crossovers, census, evaluation, gametes, loci.a, model, loci.s, auto.mp.set, sex.mp.set) {        # Makes eggs and sperm       # x = population
  GametePool <- list()
  for(i in 1:census){
    if(evaluation[i, 1] == "Male"){
      #cat("doing spermatogenesis\n")
      foo <- Spermatogenesis(x[[i]], gametes, loci.a, model, loci.s, crossovers, auto.mp.set)
    }else{
      #cat("doing oogenesis\n")
      foo <- Oogenesis(x[[i]], gametes, loci.a, model, loci.s, crossovers, auto.mp.set, sex.mp.set)
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