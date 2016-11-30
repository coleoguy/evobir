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
  if(is.matrix(x)){    
    # 1 = rows    
    # 2 = columns
    # 3 = whole matrix
    if(MARGIN == 1 | MARGIN == 2) y <- apply(x, MARGIN = MARGIN, FUN = vec.norm)
    if(MARGIN == 3) y <- matrix(vec.norm(unlist(x)), nrow(x), ncol(x))
  }
  return(y)
}