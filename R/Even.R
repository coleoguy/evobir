Even <- function(x){
  if(!is.numeric(x)){
    cat("Even only works on numbers")
    stop()
  }
  x %% 2 == 0
}
