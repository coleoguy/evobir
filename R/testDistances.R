testDistances <- function(tree, trait1, trait2, n=100, model="ARD", verbose=F){
  ###### HELPER FUNCTIONS ######
  results <- list()  
  map1 <- make.simmap(tree, x=trait1, nsim=n, model=model)
  map2 <- make.simmap(tree, x=trait2, nsim=n, model=model)
  if(n > 1){
    for(i in 1:n){
      results[[i]] <- getDist(map1[[i]], map2[[i]])
    }
  }
  if(n==1) results[[1]] <- getDist(map1, map2)
  
  Ff=Fr=Rf=Rr=fF=fR=rF=rR<-vector()
  for(i in 1:n){
    Ff <- c(Ff, results[[i]]$Ff)
    Fr <- c(Fr, results[[i]]$Fr)
    Rf <- c(Rf, results[[i]]$Rf)
    Rr <- c(Rr, results[[i]]$Rf)
    fF <- c(fF, results[[i]]$fF)
    fR <- c(fR, results[[i]]$fR)
    rF <- c(rF, results[[i]]$rF)
    rR <- c(rR, results[[i]]$rR)
  }  
  # we now have vectors describing the distances between
  # each type of change pairing
  # lets convert these to a format useful for 
  # analysis this will be a dataframe one column is time
  # (distance between changes) and the second is the type of
  # change this allows us to feed it right into something like
  # aov.
  f.result <- c(Ff, Fr, Rf, Rr, fF, fR, rF, rR)
  types <- c(rep("Ff",length(Ff)),
             rep("Fr",length(Fr)),
             rep("Rf",length(Rf)),
             rep("Rr",length(Rr)),
             rep("fF",length(fF)),
             rep("fR",length(fR)),
             rep("rF",length(rF)),
             rep("rR",length(rR)))
  df <- data.frame(f.result,types)
  pval <- summary(aov(f.result~types, data=df))[[1]][1,5]
  if(verbose==F){
    return(pval)
  }
  if(verbose==T){
    return(list(pval,df))
  }
}
