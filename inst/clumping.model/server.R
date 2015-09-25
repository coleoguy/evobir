library(shiny)


shinyServer(function(input, output) {
  
  RandomDist <- function(radius=10, 
                         individuals=5, 
                         iterations=100, 
                         observed=34){
    # this internal function just calculates a pairwise distance
    pw.dist <- function(x1, x2, y1, y2){
      sqrt(((x1 - x2) ^ 2) + ((y1 - y2) ^ 2))
    }
    # this function tests whether a point fall inside (or on)
    # a circle of given radius
    in.circle <- function(x, y, radius){
      (x - radius)^2 + (y - radius)^2 < radius^2
    }
    # we will store our final null dist here
    final.results <- vector()
    # this loop will repeat till iterations is reached
    for(k in 1:iterations){
      # this print statement allows us to monitor progress
      if(k/10 == round(k/10)) print(paste("iteration", k))
      # we will store our coordinates in this dataframe
      coords <- data.frame(x=0,y=0)
      # this will loop through for the number of individuals
      for(m in 1:individuals){
        # this flag will be false until we generate a 
        # valid coordinate
        inside <- F
        # we generate the coordinates in this while statement
        while(inside == F){
          # pick coordinates from uniform distribution
          t.coords <- runif(n=2, min = 0, max = radius*2)
          # now test if they are in circle
          if(in.circle(x=t.coords[1], y=t.coords[2], radius)){
            inside <- T
          }
        }
        coords[m, 1:2] <- t.coords
      }
      colnames(coords) <- c("X", "Y")
      # this calculates the cumulative sum of the number of 
      # individuals - 1 (this is equivelant to the number of 
      # pairwise measures)
      pairnum <- tail(cumsum((individuals - 1):1), n=1)
      pairmeas <- vector()
      counter <- 2
      # this loop is going to do the pairwise measures
      for(i in 1:(individuals-1)){
        for(j in counter:individuals){
          # here we do an actual pairwise calculation
          x <- pw.dist(coords[i,1],
                       coords[j,1], 
                       coords[i,2], 
                       coords[j,2])
          pairmeas <- c(pairmeas, x)
        }
        counter <- counter + 1
      }
      # store result from each iteration
      final.results[k] <- mean(pairmeas)
    }
    return(final.results)
  }
  
  pval.calc <- function(x, observed){
    lower <- sum(x > observed)
    higher <- sum(x < observed)
    if(lower < higher){
      pval <- lower/length(x)
    }else{
      pval <- higher/length(x)
    }
    return(pval)
  }  
  
  
  x <- reactive({
    RandomDist(radius = input$radius,
               individuals = input$individuals,
               iterations = input$iterations,
               observed = input$observed)
    
  })
  y <- reactive({
    pval.calc(x = x(),
              observed = input$observed)
 
    
  })  
  output$densPlot <- renderPlot({
      plot(density(x()),
           main = "Distribution of Mean Pairwise Measures", 
           xlab = "mean distance",
           ylab = "density",
           cex.main=1)
      abline(v=input$observed, col="red", lwd=2)
      mylims <- par("usr")
      text(x=mylims[1], 
           y=mean(mylims[3:4]), 
           paste("pvalue =", y(),
                 "\nMean Exp =", signif(mean(x()), digits=3),
                 "\nSE Null dist =", signif(sd(x())/sqrt(length(x())), digits=3)),
           col="red",
           pos=4,
           cex=1)
      
  })  
})

