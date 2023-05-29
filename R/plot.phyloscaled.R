# methods to plot a scaled phylogeny produced by scaleTreeRates

plot.phyloscaled <- function (tree,method="multiply",palette="RdYlGn",
                              edge.width=1,cex=1,show.tip.label = T)
{
  # Arguments
  # tree: tree of class phyloscaled
  # method: whether to multiply or color edges 
  # palette: diverging palette that can be passed to brewer.palette
  
  #check class
  if(!"phyloscaled" %in% class(tree)){
    print("tree is not of class phyloscaled")
    stop()
  }
  
  if(method == "multiply"){
    
    #scale edge.lengths by multiplying 
    tree$edge.length <- tree$edge.length * tree$scalar
    
    #plot
    plot.phylo(tree,edge.width = edge.width,cex = cex,
               show.tip.label=show.tip.label)
    
    
  } else if(method == "color") {
    
    #get states and number of states
    states <- sort(unique(tree$scalar))
    states.numeric <- as.numeric(as.factor(tree$scalar))
    n <- length(states)
    
    if(n >= 11){
      
      print("number of states exceeds maximum for brewer.pal, using viridis")
      
      palette <- "viridis"
      
    }
    
    #get number above and below
    slow <- sum(states < 1)
    fast <- sum(states > 1)
    
    if(slow >= 5 ||
       fast >= 5){
      
      print("number of states exceeds maximum for brewer.pal, using viridis")
      
      palette <- "viridis"
      
    }
    
    #generate colors
    if(palette=="RdYlGn"){
      if(1 %in% states){
        cols <- c(brewer.pal(n=2*slow + 1,name=palette)[1:slow],
                  brewer.pal(n=3,name=palette)[2],
                  brewer.pal(n=2*fast + 1,name = palette)[(fast + 1):
                                                            (2*fast + 1)]
        )
      } else {
        
        cols <- c(brewer.pal(n=2*slow + 1,name=palette)[1:slow],
                  brewer.pal(n=2*fast + 1,name = palette)[(fast + 2):
                                                            (2*fast + 1)]
        )
      }
    } else {
      
      if(1 %in% states){
        cols <- c(viridis(n=slow,begin=0,end=0.49),
                  viridis(n=1,begin=0.5,end=0.5),
                  viridis(n=fast,begin=0.51,end=1))
      } else {
        cols <- c(viridis(n=slow,begin=0,end=0.49),
                  viridis(n=fast,begin=0.51,end=1))
      }
      
    }
    
    #assign color to edges
    edge.cols <- c()
    
    for(i in 1:length(states.numeric)){
      edge.cols[i] <- cols[states.numeric[i]]
    }
    
    plot.phylo(tree,edge.color = edge.cols,edge.width = edge.width,cex = cex,
               show.tip.label=show.tip.label)
  } else {
    stop("unrecognized method")
  }
}
