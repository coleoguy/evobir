RadDose <- function(evaluation, sex.bias.mut) {   # list with elements of sex and index of individual in population
  x <- as.data.frame(evaluation)
  bar <- list()
  if(sex.bias.mut > sample(1:100,1)/100){
    bar[1] <- "Male"
    bar[2] <- sample(which(evaluation[,1]=="Male"), 1)
  }else{
    bar[1] <- "Female"
    bar[2] <- sample(which(evaluation[,1]!="Male"), 1)
  }
  return(bar)                                    # here we just return a list with the sex and the index number
}