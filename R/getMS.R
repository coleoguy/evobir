getMS <- function(tree, samps, report, n.site=NULL){
  # add a demes tracking column to the edge matrix
  tree$edge <- cbind(tree$edge, rep(NA, nrow(tree$edge)))
  output <- vector()
  used.demes <- 0
  check <- T
  while(check == T){
    # find tips
    tips <- tree$edge[!tree$edge[, 2] %in% tree$edge[, 1], 2]
    
    # determine which rows of #edge and elements of #edge.length these are
    branches <- which(tree$edge[,2] %in% tips)
    
    # find the shortest of the branches
    shortb <- branches[which.min(tree$edge.length[branches])]
    
    # find the parent node
    par.node <- tree$edge[shortb, 1]
    
    # find branches with this parent
    hits <- tree$edge[, 1] == par.node
    
    # remove one already stored
    hits[shortb] <- F
    
    # store the branch we will join
    shortb[2] <- which(hits)
    
    # determine if demes are named branch 1
    if(is.na(tree$edge[shortb[1], 3])){
      demes1 <- max(used.demes) + 1
      used.demes <- c(used.demes, (max(used.demes) + 1))
    }else{
      demes1 <- tree$edge[shortb[1], 3]
    }
    
    # determine if demes are named branch 1
    if(is.na(tree$edge[shortb[2], 3])){
      demes2 <- max(used.demes) + 1
      used.demes <- c(used.demes, (max(used.demes) + 1))
    }else{
      demes2 <- tree$edge[shortb[2], 3]
    }
    
    # store deme name in tree$edge
    tree$edge[tree$edge[, 2] == par.node, 3] <- demes2
    
    # record ej statement
    output <- paste(output, "-ej", tree$edge.length[shortb[1]], demes1, demes2)
    
    # check and see if we still have more work to do
    if(length(tree$edge.length)>2){
      
      #add edge length that will be lost when we remove branches
      tree$edge.length[tree$edge[, 2] == par.node] <- 
        tree$edge.length[tree$edge[, 2] == par.node] + 
        tree$edge.length[shortb[1]]
      
      
      #update tree$edge.length
      tree$edge.length <- tree$edge.length[-shortb]
      
      #update tree$edge
      tree$edge <- tree$edge[-shortb, ]
    }else{
      
      # when we are down to 2 edges we are done
      check <- F
    }
  }
  z <- paste(rep("1", length(tree$tip)), collapse=" ")
  outz <- paste("-", report, sep="")
  
  
  # added migration code paste in before "output"
  
  output <- paste(length(tree$tip), 
                  samps, "-I",
                  length(tree$tip), 
                  z, 
                  output,
                  outz,
                  n.site)
  return(output)
}
