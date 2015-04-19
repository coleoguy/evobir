Even <- function(x){
  result <- vector()
  #you should test whether its numeric and return an intelligable error
  for(i in 1:length(x)){
    if(x[i]/2 == round(x[i]/2)){
      result[i] <- T
    }else{
      result[i] <- F
    }
  }
  return(result)
}
