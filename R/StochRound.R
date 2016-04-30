StochRound <- function(x){
  ## extract the decimal portion
  q <- abs(x - trunc(x))
  
  ## draw a value 0 or 1 with probability
  ## based on how close we already are
  adj <- sample(0:1, size = 1, prob = c(1 - q, q))
  
  ## make it negative if x is
  if(x < 0) adj <- adj * -1
  
  ## return our new value
  trunc(x) + adj
}
