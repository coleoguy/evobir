library(shiny)
library(ape)
library(phytools)
shinyServer(function(input, output){

  ## Reads in the user specified tree
  getTree <- reactive({
    infile <- input$tree$datapath
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    read.tree(infile)
  })

  ## Reads in the user specified csv file  
  getData <- reactive({
    infile2 <- input$traits$datapath
    if (is.null(infile2)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    x <- read.csv(infile2, header=input$headers)
  })
    output$plot.smap <- renderPlot({
      data <- getData()
      tree <- getTree()
      if(!is.null(tree) & !is.null(data)){
        trait1 <- as.factor(data[,2])
        names(trait1) <- data[,1]
        ## stochastic mapping
        if(input$plotType == "stmap"){
          cols<-setNames(c(input$col1, input$col2),levels(trait1))
          hist <- make.simmap(tree, trait1, model = input$modelSIMMAP)
          plotSimmap(hist, 
                     colors=cols, 
                     lwd=input$lwd, 
                     fsize=input$cex.font,
                     type=input$type)
        }
        ## tip mapping
        if(input$plotType == "tips"){
          cols<-setNames(c(input$col1, input$col2),levels(trait1))
          plot(tree,
               type=input$type,
               edge.width=input$lwd,
               cex=input$cex.font)
          tiplabels(pch=input$pch, , cex=input$cex.tip, col=cols[trait1])
        }
        ## ASR mapping
        if(input$plotType == "mlasr"){
          cols<-setNames(c(input$col1, input$col2),levels(trait1))
          plot(tree,
               type=input$type,
               edge.width=input$lwd,
               cex=input$cex.font)
          hist <- ace(trait1, tree, type="discrete",
                      model = input$modelACE)
          nodelabels(pie=hist$lik.anc[,1], 
                     cex=input$cex.tipML, piecol=c(input$col1, input$col2))
        }
      }
    })
  })
