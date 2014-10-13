library(shiny)
shinyServer(function(input, output) {
  tree <- reactive({
    x <- seq(input$N.mean - 4 * input$N.sd, 
             input$N.mean + 4 * input$N.sd, 
             length = 200)
    y <- dnorm(x, mean = input$N.mean, sd=input$N.sd)
  })
  output$treePlot <- renderPlot({
      plot(x=tree()$x, y=tree()$y, col = "red", type="l", lwd=3)
  })  
})

