library(shiny)


shinyServer(function(input, output) {
  
  x <- reactive({
    if(input$select == 1){
      seq(input$mu - 4 * input$sigma, 
          input$mu + 4 * input$sigma, 
          length = 200)
    }else if(input$select == 2){
      high <- qexp(.999, rate = input$lambda, lower.tail = TRUE, log.p = FALSE)
      seq(0, high, length = 200)
    }else if(input$select == 3){
      high <- qgamma(.999, shape = input$kappa, rate = input$theta)
      seq(0, high, length = 200)
    }else if(input$select == 4){
      low <- qlogis(.01, location = input$mu, scale = input$sigma*0.551328895, lower.tail = TRUE, log.p = FALSE)
      high <- qlogis(.99, location = input$mu, scale = input$sigma*0.551328895, lower.tail = TRUE, log.p = FALSE)
      seq(low, high, length = 200)      
    }
  })
  y <- reactive({
    if(input$select == 1){
      dnorm(x(), mean = input$mu, sd=input$sigma)
    }else if(input$select == 2){
      dexp(x(), rate = input$lambda)
    }else if(input$select == 3){
      dgamma(x(), shape = input$kappa, rate = input$theta)
    }else if(input$select == 4){
      dlogis(x(), location = input$mu, scale = input$sigma*0.551328895)
  }
})
   
  output$treePlot <- renderPlot({
      plot(x=x(), y=y(), col = "red", ylab="density", 
           xlab=paste("Distribution:"), type="l", lwd=3)
  })  
})

