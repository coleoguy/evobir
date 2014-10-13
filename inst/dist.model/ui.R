# Define UI for bd trees
shinyUI(pageWithSidebar(
  headerPanel("Probability Density of Distributions"),    # Application title
  sidebarPanel(
    selectInput("select", label = "Select Statistical Distribution",
                choices = list("Normal" = 1, "Exponential" = 2,
                               "Gamma" = 3, "Logistic" = 4), selected = 1), 
    sliderInput("mu", "mu:", 
                min = -100, max = 100, value=0, step = 5),
    sliderInput("sigma", "sigma:", 
                min = 0.1, max = 10, value=1, step = .1),
    sliderInput("lambda", "lambda:", 
                min = 0.1, max = 2, value=1, step = .01),
    sliderInput("kappa", "kappa:", 
                min = 0.1, max = 10, value=1, step = .1),
    sliderInput("theta", "theta:", 
                min = 0.3, max = 10, value=1, step = .1)
  ),  
  mainPanel(
    plotOutput("treePlot", width='100%', height='600px')
  )
))