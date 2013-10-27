library(shiny)


shinyServer(function(input, output) {
  tree <- reactive({
   ## put the call to create the tree here
  })
  

  output$genePlot <- renderPlot({
   ## put the call to create the plot here
    
  })
  
})


