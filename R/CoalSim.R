CoalSim <- function(census = 15, lw = 2, ln.col = 'blue'){   
    history <- as.data.frame(1:census)
    lineages <- census
    counter <- 1
    counter <- counter + 1
    history[1:census, counter] <- sample(1:census, length(unique(history[, (counter - 1)])), replace = T)
    while(lineages > 1){
      lineages <- sum(!is.na(unique(history[1:10, counter])))
      history[1:lineages, (counter + 1)] <- sort(unique(history[1:10, counter]))
      counter <- counter + 1
      history[1:lineages, (counter + 1)] <- sort(sample(1:census, length(unique(history[, (counter)])) - 1, replace = T))
      counter <- counter + 1
    }
    ## now plot it
    plot(0, 0, col = 'white', ylim = c(1,census),
         xlim = c(1, counter/2 + 1),
         main='Coalesence Simulation', 
         xlab = 'Generation', 
         ylab = 'Individuals')
    foo <- 1
    foox <- 1
    while(foo < (counter)){
      for(i in 1:census){  
        lines(c(foox,(foox + 1)), history[i,foo:(foo + 1)], lwd = lw, col=ln.col)
      }
      foo <- foo + 2
      foox <- foox + 1
    }
}