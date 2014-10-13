library(shiny)
shinyServer(function(input, output) {
  tree <- reactive({
    x <- seq(input$N.mean - 4 * input$N.sd, 
             input$N.mean + 4 * input$N.sd, 
             length = 200)
    y <- dnorm(x, mean = N.mean, N.sd)
  })
  counts <- 
  output$treePlot <- renderPlot({
      plot(tree())
  })  
})


