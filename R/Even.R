Even <- function(x){
  result <- vector()
  for(i in 1:length(x)){
    if(x[i]/2 == round(x[i]/2)){
      result[i] <- T
    }else{
      result[i] <- F
    }
  }
  return(result)
}