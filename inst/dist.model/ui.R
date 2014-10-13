# Define UI for bd trees
shinyUI(pageWithSidebar(
  headerPanel("Probability Density of Distributions"),    # Application title
  sidebarPanel(
    tags$div(class="header", checked=NA,
             tags$p("Normal Distribution")),
    tags$div(class="header", checked=NA,
             tags$p("Select Statistical Distribution")),
    selectInput("select", label = h3("Select box"),
                choices = list("Normal" = 1, "Exponential" = 2,
                               "Gamma" = 3, "Logistic" = 4), selected = 1),    
    
    sliderInput("time", "Time for simulation:", 
                min = 4, max = 20, value=10, step = 1),
    sliderInput("birth", "Birth Rate:", 
                min = .125, max = .75, value=.5, step = .125),
    sliderInput("death", "Death Rate:", 
                min = 0, max = .25, value=0, step = .05),
    checkboxInput(inputId = "extinct",
                  label = "Prune Extinct Taxa",
                  value = FALSE),
    actionButton("seed.val", 'Simulate')
    
  ),  
  mainPanel(
    plotOutput("treePlot", width='100%', height='600px')
  )
))