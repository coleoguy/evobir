countTrees <- function(collection = NULL, ref = NULL, verbose=T){
  if(is.null(collection)) stop("supply a path to a collection of trees in a Newick format file")
  if(is.null(ref)) stop("supply a path to a set of topologies to count in a Newick format file")
  trees <- read.tree(collection)
  types <- read.tree(ref)
  top.num <- length(types)
  class(types) <- "multiPhylo"
  classification <- vector(length=top.num)
  bad <- c()
  for(i in 1:length(trees)){
    missing <- T
    counter <- 1
    while(missing){
      x <- all.equal.phylo(trees[[i]], types[[counter]], use.edge.length=F)
      if(x){
        classification[counter] <- classification[counter] + 1
        missing <- F
      }else if(counter == top.num){
        bad <- c(bad, i)
        missing <- F
      }
      counter <- counter + 1
    }
  }
  if(verbose){
    if(length(bad)>0){
      print(paste("Some trees do match available topologies. You may want to check trees:", bad))
    }
  }
  return(classification)
}
