library(shiny)


shinyServer(function(input, output) {
  
  x <- reactive({
    if(input$select == 1){
      seq(input$mu - 40, 
          input$mu + 40, 
          length = 200)
    }else if(input$select == 2){
      high <- qexp(.999, rate = input$lambda, lower.tail = TRUE, log.p = FALSE)
      seq(0, high, length = 200)
    }else if(input$select == 3){
      high <- qgamma(.999, shape = input$kappa, rate = input$theta)
      seq(0, high, length = 200)
    }else if(input$select == 4){
      seq(input$mu - 30, 
          input$mu + 30, 
          length = 200)
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
    if(input$select == 1){
      plot(x=x(), y=y(), col = "red", ylab="density", ylim=c(0,.6),
           xlab=paste("Distribution:"), type="l", lwd=3)
      abline(h=0, lty=3, cex=2)
    }else if(input$select == 2){
      plot(x=x(), y=y(), col = "red", ylab="density", ylim=c(0,2),
           xlab=paste("Distribution:"), type="l", lwd=3)
      abline(h=0, lty=3, cex=2)
    }else if(input$select == 3){
      plot(x=x(), y=y(), col = "red", ylab="density", 
           xlab=paste("Distribution:"), type="l", lwd=3)
      abline(h=0, lty=3, cex=2)
    }else if(input$select == 4){
      plot(x=x(), y=y(), col = "red", ylab="density", 
           xlab=paste("Distribution:"), type="l", lwd=3)
      abline(h=0, lty=3, cex=2)
    }
  })  
})

