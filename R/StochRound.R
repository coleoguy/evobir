StochRound <- function(x){
  ## extract the decimal portion
  q <- abs(x - trunc(x))
  ## draw a value 0 or 1 with probability
  ## based on how close we already are
  adj <- c()
  for(i in 1:length(x)){
    adj[i] <- sample(0:1, size = 1, prob = c(1 - q[i], q[i]))
    ## make it negative if x is
    if(x[i] < 0) adj[i] <- adj[i] * -1
  }
  ## return our new value
  trunc(x) + adj
}
