library(shiny)
library(phytools)
shinyServer(function(input, output) {
  tree <- reactive({
      set.seed(input$seed.val)
      pbtree(b = input$birth, 
             d = input$death,
             t = input$time, 
             scale = 1,
             nsim = 1,
             extant.only = input$extinct) 
  })
  counts <- 
  output$treePlot <- renderPlot({
      plot.phylo(tree()[[i]], show.tip.label=F)
      mtext(paste("N =", length(tree()[[i]]$tip.label)), side = 1, line = 0)
  })  
})


