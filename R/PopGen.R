
PopGen <- function(fitness = c(.98, 2, .98), initial.A = .5, pop = 1000, 
                   gen = 100, var.plot = 1, iter = 20,
                   chrom = "Auto", width){
  genotypes <- c('AA', 'Aa', 'aa')
  iter.col <- rainbow(iter)
  plot(0,0, xlab = "Generation", xlim = c(0, gen),
       ylab = paste('Frequency of', genotypes[var.plot], 'genotype'),
       col='white', ylim = c(0,1),
       main= paste('Fitness of', genotypes[var.plot], '=', fitness[var.plot] ))
  for(k in 1:iter){
    adults <- c(rep(1, each = round(pop*initial.A^2)), 
                rep(2, each = round(pop*2*initial.A*(1-initial.A))), 
                rep(3, each = round(pop*(1-initial.A)^2)))
    plot.val <- vector()
    for(i in 1:gen){
      A <- ((2 * sum(adults == 1)) + sum(adults ==2) ) / (pop*2)
      babies <-  c(rep(1, each = round(pop*A^2)), 
                   rep(2, each = round(pop*2*A*(1-A))), 
                   rep(3, each = round(pop*(1-A)^2)))
      pop.fit <- vector() #fitness for each offspring
      for(j in 1:length(babies)){
        if(babies[j] == 1){
          pop.fit[j] <- fitness[1]
        }else if(babies[j] == 2){
          pop.fit[j] <- fitness[2]
        }else{
          pop.fit[j] <- fitness[3]
        }
        }
      adults <- sample(babies, pop, replace = T, prob = pop.fit)
      if(var.plot == 1){
        plot.val[i] <- sum(adults == 1)
      }else if(var.plot == 2){
        plot.val[i] <- sum(adults == 2)
      }else{
        plot.val[i] <- sum(adults == 3)
      }
    }
    lines(plot.val/pop, lwd = width, col=iter.col[k])
  }
}
  
