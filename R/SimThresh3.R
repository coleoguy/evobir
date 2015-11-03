SimThresh3 <- function(tree, liabilities=F){
  returnState <- function(x,y){
    state <- vector()
    # angle A is at origin
    # angle B is at 1,0
    # angle C is at x,y
    # so we can calculate the sides like this
    a <- sqrt((x-1)^2 + (y)^2)
    c <- sqrt(x^2+y^2)
    b <- sqrt(1)
    # and then we can calculate the angle at the origin like this:
    A <- (acos((b^2 + c^2 - a^2)/(2*b*c))) * (180/pi)
    state <- vector(mode="numeric", length=length(y))
    state[(y >= 0 & A < 90) | (y < 0 & A < 30)] <- 1
    state[(y >= 0 & A > 90) | (y < 0 & A > 150)] <- 2
    state[y < 0 & (A > 30 & A < 150)] <- 3
    names(state) <- names(x)
    return(state)
  }
  x<-fastBM(tree)
  y<-fastBM(tree)
  tip.state <- returnState(x,y)
  if(liabilities == F) return(tip.state)
  if(liabilities == T){
    results <- vector("list", 3)
    results[[1]] <- tip.state
    results[[2]] <- x
    results[[3]] <- y
    names(results) <- c("observed", "liab1", "liab2")
    return(results)
  } 
}
