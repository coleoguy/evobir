library(shiny)
shinyServer(function(input, output) {
  x <- reactive({
    if(select == 1){
      seq(input$N.mean - 4 * input$N.sd, 
          input$N.mean + 4 * input$N.sd, 
          length = 200)
    }
    if(select == 2){
      seq(input$N.mean - 6 * input$N.sd, 
          input$N.mean + 6 * input$N.sd, 
          length = 200)
    }
    
    
    
  })
  y <- reactive({
    dnorm(x(), mean = input$N.mean, sd=input$N.sd)
  })
  output$treePlot <- renderPlot({
      plot(x=x(), y=y(), col = "red", type="l", lwd=3)
  })  
})

