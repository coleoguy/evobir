library(shiny)
shinyServer(function(input, output) {
 # tree <- reactive({
    #if(select == 1){
    #  x <- seq(N.mean - 4 * N.sd, N.mean + 4 * N.sd, length = 200)
    #  y <- dnorm(x, mean = N.mean, N.sd)
    #}
      
 # })
  output$treePlot <- renderPlot({
    plot(1:10,1:10)
#    if(select == 1){
#      tree(N.mean, N.sd)
#    }
#    plot(x=x,y=y,col="red",  type="l", ylab="Density")

    
  })  
})


