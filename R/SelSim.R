SelSim <- function(census, initial.freq, selection.coef, iter, generations){
  results <- matrix(,iter, generations)
  for(j in 1:iter){
    population <- vector()
    freq <- vector()
    count1 <- round(census * initial.freq)
    count0 <- census - count1
    population[1:census] <- c(rep(0, count0), rep(1, count1))
    for(i in 1:generations){
      pop.prob<- vector()
      for(k in 1:census){
        if(population[k] == 0){
          pop.prob[k] <- 1
        }else{
          pop.prob[k] <- selection.coef
        }
      }
      population[1:census] <- sample(population, census, prob = pop.prob, replace = T)
      freq[i] <- sum(population == 1) / census
    }
    results[j, ] <- freq
  }
  plot(0, 0, ylim = c(0,1), xlim = c(0,generations), col = 'white', 
       main = paste("Change in Freq of a Mutation\n initial = ", initial.freq,
                    "  S = ", selection.coef, "  N = ", census, sep = "" ),
       xlab = 'Generations',
       ylab = 'Allele Frequency')
  gene.cols <- rainbow(iter)
  for(i in 1:iter) lines(results[i, 1:generations], col = gene.cols[i], lwd = 1.5)
  fixed <- sum(results[,generations] == 1)
  lost <- sum(results[,generations] == 0)
  segregating <- sum(results[,generations] > 0 & results[,generations] < 1)
  if(fixed > 0){
    if(fixed > 1) print(paste(fixed, 'Iterations fixed the selected mutation.'))
    if(fixed == 1) print(paste(fixed, 'Iteration fixed the selected mutation.'))
  }
  if(lost > 0){
    if(lost > 1) print(paste(lost, 'Iterations lost the selected mutation.'))
    if(lost == 1) print(paste(lost, 'Iteration lost the selected mutation.'))
  }
  if(segregating > 0){
    if(segregating > 1) print(paste(segregating, 'Iterations still have the selected mutation segregating.'))
    if(segregating == 1) print(paste(segregating, 'Iteration still have the selected mutation segregating.'))
  }
}