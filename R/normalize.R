normalize <- function(x, MARGIN=1){
  vec.norm <- function(x){
    x[1:length(x)] <- ((x - min(x)) / (max(x) - min(x)))
    return(x)
  }
  if(is.vector(x, mode="numeric")){
    y <- vec.norm(x)
  }
  if(is.vector(x, mode="list")){
    y <- lapply(x, vec.norm)
  }
  if(is.matrix(x)){    # 1 = rows    2 = columns
    y <- apply(x, MARGIN = MARGIN, FUN = vec.norm)
  }
  return(y)
}