library(shiny)
# I like this color pallete much more than the base options
library(viridis)
shinyServer(function(input, output) {
  # this is the actual brownian motion simulation
  x <- reactive({
    # insure reproducability
    set.seed <- input$seed.val
    replicate(input$reps, cumsum(rnorm(n=input$gens, sd=input$rate)))
  })
  # Now we will build the output plot
  output$distPlot <- renderPlot({  
    # we will have two plots
    par(mfcol=c(2,1))
    # need a bit more room 
    par(mar=c(3,2,1,1))
    # this just keeps our color variation from being hidden when itterations++
    cols <- sample(viridis(input$reps))
    # the primary plot
    matplot(x(),type="l",lty=1, lwd=1.5, col=cols, main="Brownian Motion", ylab="Trait Value", xlab="Generations")
    # this line will indicate what is being plotted in the histogram below
    abline(v=input$mon, lwd=3, col="red")
    # now we plot the distribution switched from a histrogram this looks prettier IMHO
    plot(density(x()[input$mon,]), 
         xlim=c(min(x()),max(x())),
         col = 'red', 
         xlab="Trait value",
         ylab="",
         lwd=3,
         main=paste("Trait distribution at generation", input$mon))
    
  })
})
