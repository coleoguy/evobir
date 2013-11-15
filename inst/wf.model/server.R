## right now I am not getting fixations of A if population size fluctuates?


library(shiny)
genotype.options <- c('AA', 'Aa', 'aa', 'A', 'a')
gen.exp <- function(x, y, wAA, wAa, waa, qAa, qaA){
  foo <- vector()
  A <- x
  AA <- A^2       #p2
  Aa <- 2*A*(1-A) #2pq
  aa <- (1-A)^2   #q2
  for(i in 1:y){
    w.bar <- AA*wAA + Aa*wAa + aa*waa  #mean fitness
    AA <- AA * (wAA / w.bar)  
    Aa <- Aa * (wAa / w.bar)
    aa <- aa * (waa / w.bar)
    foo[i] <- A <- (AA + .5*Aa)
    if(qAa + qaA != 0) A <- A + {{1 - A} * qaA} - {A * qAa} # 
    AA <- A^2
    Aa <- 2*A*(1-A)
    aa <- (1-A)^2
  }
  return(foo)
}
ShinyPopGen <- function(fitness, initial.A, pop, gen, var.plot, iter, heath, qAa, qaA, popD){
  results <- matrix(,iter,gen)
  pop2 <- vector()
  for(k in 1:iter){                             # this loop goes through the iterations
    adults <- c(rep(1, each = round(pop*initial.A^2)), 
                rep(2, each = round(pop*2*initial.A*{1-initial.A})), 
                rep(3, each = round(pop*{1-initial.A}^2)))
    plot.val <- vector()
    for(i in 1:gen){                            # this loop goes through the generations
      if((2 * sum(adults == 1) + sum(adults ==2)) / {pop*2} < 1){
        A <- (2 * sum(adults == 1) + sum(adults ==2)) / {pop*2}
      }else{
        A <- 1
      }
      if(qAa + qaA != 0) A <- A + {{1 - A} * qaA} - {A * qAa}
      pop2 <- pop
      if(popD != 0) pop2 <- pop + runif(1, min = -popD, max= popD)
      if(pop2 < 10) pop2 <- 10
      if(A>1) A <- 1
      #print(paste("gen:",i, "\npop2:",pop2, "A:", A))
      babies <-  c(rep(1, each = round(pop2*A^2)), 
                   rep(2, each = round(pop2*2*A*{1-A})), 
                   rep(3, each = round(pop2*(1-A)^2)))
      pop.fit <- vector(length = length(babies))                       # fitness for each offspring
      pop.fit[babies == 1] <- fitness[1]
      pop.fit[babies == 2] <- fitness[2]
      pop.fit[babies == 3] <- fitness[3]
      adults <- sample(babies, pop2, replace = T, prob = pop.fit)
      AA <- sum(adults == 1)
      Aa <- sum(adults == 2)
      plot.val[i] <- AA + .5 * Aa
    }
    results[k,]<-plot.val
  }
  return(results)
}
fates <- function(data.table, generations, population, iterations){
  lost <- 0
  fixed <- 0
  for(i in 1:iterations){
    if(data.table[i,generations] == (2 * population)) fixed <- fixed + 1
    if(data.table[i,generations] == 0) lost <- lost + 1
  }  
  return(c(lost, fixed))
}
shinyServer(function(input, output) {
  genotypes <- reactive({
    paste('Frequency of', genotype.options[as.numeric(input$var.plot)])
    })
  
  data <- reactive({   # this will contain the number of A alleles in the population
    set.seed <- input$seed.val
    ShinyPopGen(fitness=c(input$fit.AA, input$fit.Aa, input$fit.aa), 
                                initial.A = input$initial.A, 
                                pop = input$pop, 
                                gen = input$gen, 
                                var.plot = input$var.plot, 
                                iter = input$iter,
                                heath = input$heath, input$qAa, input$qaA,
                                input$popD)
  })
  
  fate.lost <- reactive({sum(data()[,input$gen] == 0)})
  fate.fixed <- reactive({sum(data()[,input$gen] == (input$pop))})
  
  
  output$caption1 <- renderText({paste('Allele A lost in', fate.lost(), 'populations')})
  output$caption2 <- renderText({paste('Allele A fixed in', fate.fixed(), 'populations')})
  
  expected.A <- reactive({
    gen.exp(input$initial.A, input$gen, input$fit.AA, input$fit.Aa, input$fit.aa,
            input$qAa, input$qaA)
  })
  
  output$genePlot <- renderPlot({
    if(input$plotpop == F){
      plot(0, 0, col = 'white', ylim = c(0, 1), xlim = c(0, input$gen),
         xlab = 'Generations', ylab = genotypes(), cex.lab=1.2)
      mtext("Produced with the package evobiR", side = 1, cex=.8, line=4)
      for(i in 1:input$iter){
        if(input$var.plot == 1){
          lines(1:input$gen, (data()[i,1:input$gen]/input$pop)^2,
          col=rainbow(input$iter)[i], lwd = input$width)
        }else if(input$var.plot == 2){
          lines(1:input$gen, 2 * (data()[i,1:input$gen]/input$pop) * (1-(data()[i,1:input$gen]/input$pop)),
          col=rainbow(input$iter)[i], lwd = input$width)
        }else if(input$var.plot == 3){
          lines(1:input$gen, (1-(data()[i,1:input$gen]/input$pop))^2,
          col=rainbow(input$iter)[i], lwd = input$width)
        }else if(input$var.plot == 4){
          lines(1:input$gen, data()[i,1:input$gen]/input$pop,
          col=rainbow(input$iter)[i], lwd = input$width)
        }else{
          lines(1:input$gen, 1 - data()[i,1:input$gen]/input$pop,
          col=rainbow(input$iter)[i], lwd = input$width)
        }
      }
      if(input$traj == T){
        if(input$var.plot == 1) lines(1:input$gen, expected.A()^2, col='black', lwd = input$width+1, lty=2)
        if(input$var.plot == 2) lines(1:input$gen, 2 * expected.A() * (1-expected.A()), col='black', lwd = input$width+1, lty=2)
        if(input$var.plot == 3) lines(1:input$gen, (1-expected.A())^2, col='black', lwd = input$width+1, lty=2)
        if(input$var.plot == 4) lines(1:input$gen, expected.A(), col='black', lwd = input$width+1, lty=2)
        if(input$var.plot == 5) lines(1:input$gen, (1-expected.A()), col='black', lwd = input$width+1, lty=2)
      }
    }else{
      plot(0, 0, col = 'white', ylim = c(-.5, 1), xlim = c(0, input$gen),
           xlab = 'Generations', yaxt="n", ylab = genotypes(), cex.lab=1.2)
      mtext("Produced with the package evobiR", side = 1, cex=.8, line=4)
      for(i in 1:input$iter){
        if(input$var.plot == 1){
          lines(1:input$gen, (data()[i,1:input$gen]/input$pop)^2,
                col=rainbow(input$iter)[i], lwd = input$width)
        }else if(input$var.plot == 2){
          lines(1:input$gen, 2 * (data()[i,1:input$gen]/input$pop) * (1-(data()[i,1:input$gen]/input$pop)),
                col=rainbow(input$iter)[i], lwd = input$width)
        }else if(input$var.plot == 3){
          lines(1:input$gen, (1-(data()[i,1:input$gen]/input$pop))^2,
                col=rainbow(input$iter)[i], lwd = input$width)
        }else if(input$var.plot == 4){
          lines(1:input$gen, data()[i,1:input$gen]/input$pop,
                col=rainbow(input$iter)[i], lwd = input$width)
        }else{
          lines(1:input$gen, 1 - data()[i,1:input$gen]/input$pop,
                col=rainbow(input$iter)[i], lwd = input$width)
        }
      }
      if(input$traj == T){
        if(input$var.plot == 1) lines(1:input$gen, expected.A()^2, col='black', lwd = input$width+1, lty=2)
        if(input$var.plot == 2) lines(1:input$gen, 2 * expected.A() * (1-expected.A()), col='black', lwd = input$width+1, lty=2)
        if(input$var.plot == 3) lines(1:input$gen, (1-expected.A())^2, col='black', lwd = input$width+1, lty=2)
        if(input$var.plot == 4) lines(1:input$gen, expected.A(), col='black', lwd = input$width+1, lty=2)
        if(input$var.plot == 5) lines(1:input$gen, (1-expected.A()), col='black', lwd = input$width+1, lty=2)
      }
      abline(h=0)
    text(1,-.1, "I am working on this", col = "red",cex=2, pos=4)
    #axis(side = 4)
    #mtext(side = 4, line = 3, "Population Size")
    }
})

  })



