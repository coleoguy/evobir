getPaths <- function(tree){
  # extract the table
  tab <- tree$edge
  # id the tips
  tips <- tab[,2][!tab[,2] %in% tab[,1]]
  # id the root
  root <- unique(tab[,1][!tab[,1] %in% tab[,2]])
  # create a list to store results in
  paths <- list()
  # a loop to go through all tips
  for(i in 1:length(tips)){
    # this is the little dynamic portion
    x <- tips[i]
    while(x[length(x)] != root){
      y <- tab[which(x[length(x)] == tab[,2]), 1]
      x[(length(x)+1)] <- y
    }
    # lets reverse it because I think root to tip
    paths[[i]] <- rev(x)
  }
  return(paths)
}
