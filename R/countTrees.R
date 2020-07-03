countTrees <- function(collection = NULL, ref = NULL, classes=T, verbose=T){
  if(is.null(collection)) stop("supply a path to a collection of trees in a Newick format file")
  if(is.null(ref)) stop("supply a path to a set of topologies to count in a Newick format file")
  trees <- read.tree(collection)
  types <- read.tree(ref)
  top.num <- length(types)
  tree.class <- vector(length=length(trees))
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
        tree.class[i] <- counter
        missing <- F
      }else if(counter == top.num){
        bad <- c(bad, i)
        tree.class[i] <- NA
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
  if(classes==T){
    classification <- list(classification, tree.class)
    return(classfication)
  }else{
    return(classification)
  }
}
