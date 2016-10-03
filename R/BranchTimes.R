BranchTimes <- function(x){
  bt <- branching.times(x)
  n.bt <- vector(length=nrow(x$edge))
  for(i in 1:nrow(x$edge)){
    n.bt[which(x$edge[,1] == Ntip(x)+i)] <- bt[i]
  }
return(n.bt)
}
