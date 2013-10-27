library(shiny)
shinyServer(function(input, output) {
  tree <- reactive({
   ## put the call to create the tree here
    pbtree(b = input$birth, 
           n = input$num.taxa, 
           scale = 1,
           extant.only = input$extinct) 
  })
  output$genePlot <- renderPlot({
    plot.phylo(tree, show.tip.label=F)
  })  
})


