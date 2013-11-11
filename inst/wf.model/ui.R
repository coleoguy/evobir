# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Change in Genotype Frequency"),
  
  sidebarPanel(
    sliderInput("initial.A", "Initial Frequency of A allele:", 
                min = 0, max = 1, value=.5, step =.025),
    selectInput("pop", "Population Size:", 
                list(10=10,
                     50=50,
                     100=100,
                     200=200,
                     400=400))
                     min = 10, max = 500, value = 100, step = 1),
    sliderInput("gen", "Generations to simualte:", 
                min = 10, max = 500, value = 100, step = 1),
    sliderInput("fit.AA", "Fitness of AA:", 
                min = 0, max = 2, value = 1, step = .05),
    sliderInput("fit.Aa", "Fitness of Aa:", 
                min = 0, max = 2, value = 1, step = .05),
    sliderInput("fit.aa", "Fitness of aa:", 
                min = 0, max = 2, value = 1, step = .05),
    sliderInput("iter", "Iterations:", 
                min = 1, max = 50, value = 10, step = 1),
    sliderInput("width", "Line width:", 
                min = .2, max = 6, value = 2, step = .1),
    selectInput("var.plot", "Plot:",
                list("A" = 4, 
                     "a" = 5, 
                     "AA" = 1,
                     "Aa"  = 2,
                     "aa"  = 3)),
    selectInput("heath", "Benchmarking Tests:",
                list("option 1" = "preset", 
                     "option 2" = "fly")),
    checkboxInput(inputId = "traj",
                  label = "Show expected outcome",
                  value = FALSE),
    actionButton("seed.val", 'Refresh')
  ),  
  mainPanel(
       h4(textOutput("caption1")),
       h4(textOutput("caption2")),
       plotOutput("genePlot", width='100%', height='600px')
  )
))