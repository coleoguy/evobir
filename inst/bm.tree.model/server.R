library(shiny)
# I like this color pallete much more than the base options
library(viridis)

history <- function(gens, birth){
  origins <- 1
  lineages <- list()
  # the root of the tree
  lineages[[1]] <- 0
  # loop through generations with i
  for(i in 2:(gens)){
    # create next gen of every existing lineage
    for(j in 1:length(lineages)){
      lineages[[j]] <- c(lineages[[j]], rnorm(1))
    }
    # give each lineage an opportunity to split/speciate
    for(k in 1:length(lineages)){
      if(runif(1, min=0, max=1) < birth){
        lineages[[length(lineages) + 1]] <- cumsum(lineages[[k]])[length(lineages[[k]])]
        origins <- c(origins, i)
      }
    }
  }
  results <- list()
  results[[1]] <- lineages
  results[[2]] <- origins
  return(results)
}

shinyServer(function(input, output) {
  # this is the actual brownian motion simulation
  x <- reactive({
    # insure reproducability
    set.seed <- input$seed.val
    history(input$gens, input$birth)
  })
  # Now we will build the output plot
  output$distPlot <- renderPlot({  
    par(mfcol=c(1,2))
    lineages <- x()[[1]]
    origins <- x()[[2]]
    max <- min <- vector()
    for(j in 1:length(lineages)){
      max[j] <- max(cumsum(lineages[[j]]))
      min[j] <- min(cumsum(lineages[[j]]))
    }
    max <- max(max)
    min <- min(min)
    plot(x=1:(input$gens),y=cumsum(lineages[[1]]),ylim=c(min,max), type="l",
         main="BM - birth process", xlab="generations",
         ylab="trait values")
    for(j in 1:length(lineages)){
      lines(x=origins[j]:input$gens,y=cumsum(lineages[[j]]), col=viridis(length(lineages))[sample(length(lineages),1)])
    }
    abline(v=input$mon, lwd=3, col="red")
    time.slice <- vector()
    for(i in 1:length(lineages)){
      if(origins[i]<input$mon){
        time.slice <- c(time.slice, cumsum(lineages[[i]])[input$mon-origins[i]+1])
      }
    }
    time.slice <- time.slice[!is.na(time.slice)]
    # now we plot the distribution switched from a histrogram this looks prettier IMHO
    hist(time.slice, 
         col = 'red', 
         xlab="Trait value",
         ylab="",
         main=paste("Trait distribution at generation", input$mon))
  })
})
