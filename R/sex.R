sex <- function(x) {   # this function checks the sex of a genome
  if(x[7, 1] == 0){ # we dont check the Y here that way we can have XO males
    foo <- "Male"
  }else{
    foo <- "Female"
  }
  return(foo)
}
