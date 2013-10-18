DriftSim <- function(census, initial.freq, iter, generations){
  results <- matrix(,iter, generations)
  for(j in 1:iter){
    population <- vector()
    freq <- vector()
    count1 <- round(census * initial.freq)
    count0 <- census - count1
    population[1:census] <- c(rep(0, count0), rep(1, count1))
    for(i in 1:generations){
      population[1:census] <- sample(population, census, replace = T)
      freq[i] <- sum(population == 1) / census
    }
    results[j, ] <- freq
  }
  plot(0, 0, ylim = c(0,1), xlim = c(0,generations), col = 'white', 
       main = paste("Change in Freq of a Neutral Mutation\n initial freq. = ", initial.freq, 
                    "N = ", census, sep = "" ),
       xlab = 'Generations',
       ylab = 'Allele Frequency')
  gene.cols <- rainbow(iter)
  for(i in 1:iter) lines(results[i, 1:generations], col = gene.cols[i], lwd = 1.5)
  fixed <- sum(results[,generations] == 1)
  lost <- sum(results[,generations] == 0)
  segregating <- sum(results[,generations] > 0 & results[,generations] < 1)
  if(fixed > 0){
    if(fixed > 1) print(paste(fixed, 'Iterations fixed the neutral mutation.'))
    if(fixed == 1) print(paste(fixed, 'Iteration fixed the neutral mutation.'))
  }
  if(lost > 0){
    if(lost > 1) print(paste(lost, 'Iterations lost the neutral mutation.'))
    if(lost == 1) print(paste(lost, 'Iteration lost the neutral mutation.'))
  }
  if(segregating > 0){
    if(segregating > 1) print(paste(segregating, 'Iterations still have the neutral mutation segregating.'))
    if(segregating == 1) print(paste(segregating, 'Iteration still have the neutral mutation segregating.'))
  }
  top <- paste("Fixed: ", fixed)
  bot <- paste("Lost: ", lost)
  text(0, 1, top, cex=.75, pos=4)
  text(0, 0, bot, cex=.75, pos=4)
}