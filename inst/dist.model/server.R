library(shiny)
shinyServer(function(input, output) {
  tree <- reactive({
    if(select == 1){
      x <- seq(N.mean - 4 * N.sd, N.mean + 4 * N.sd, length = 200)
      y <- dnorm(x, mean = N.mean, N.sd)
    }
      
  })
  counts <- 
  output$treePlot <- renderPlot({
    tree()
    plot(x=x,y=y,col="red",  type="l", ylab="Density")

    
  })  
})


