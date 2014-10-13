# Define UI for bd trees
shinyUI(pageWithSidebar(
  headerPanel("Probability Density of Distributions"),    # Application title
  sidebarPanel(
    selectInput("select", label = "Select Statistical Distribution",
                choices = list("Normal" = 1, "Exponential" = 2,
                               "Gamma" = 3, "Logistic" = 4), selected = 1), 
    sliderInput("N.mean", "mean:", 
                min = -100, max = 100, value=0, step = 5),
    sliderInput("N.sd", "standard deviation:", 
                min = 0, max = 10, value=1, step = .25),
    sliderInput("E.lambda", "Rate parameter (lambda):", 
                min = 0, max = 5, value=1, step = .25),
    sliderInput("G.shape", "Shape parameter:", 
                min = 0, max = 10, value=1, step = .25),
    sliderInput("G.scale", "Scale parameter:", 
                min = 0, max = 10, value=1, step = .25)
  ),  
  mainPanel(
    plotOutput("treePlot", width='100%', height='600px')
  )
))