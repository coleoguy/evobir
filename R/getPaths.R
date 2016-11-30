getPaths <- function(tree, type) {
  # extract the table
  tab <- tree$edge
  # id the tips - occur in right column but not left
  tips <- tab[, 2][!tab[, 2] %in% tab[, 1]]
  
  # could look at vector of node labels
  # should be just ntips
  
  # id the root - occurs in left column but not right
  root <- (tab[, 1][!tab[, 1] %in% tab[, 2]])[1]
  # create a list to store results in
  paths <- list()
  # a loop to go through all tips
  for (i in 1:length(tips)) {
    # pull the current tip
    x <- tips[i]
    # vector to store branches in
    if(type == "branch") bl <- vector()
    # checks to see if we found root yet
    while (x[length(x)] != root) {
      # which node leads to the most recently sampled node
      y <- tab[which(x[length(x)] == tab[, 2]), 1]
      # the row of tab is equivelant to the branch index so
      # we get our branch id this way
      if(type == "branch"){
        bl <- c(bl, which(x[length(x)] == tab[, 2]))
      }
      # store our new node
      x[(length(x) + 1)] <- y
    }
    # lets reverse it because I think root to tip
    if (type == "node") paths[[i]] <- rev(x)
    if (type == "branch") paths[[i]] <- rev(bl)
  }
  return(paths)
}
