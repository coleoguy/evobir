# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Change in Genotype Frequency"),
  
  sidebarPanel(
    sliderInput("initial.A", "Initial Frequency of A allele:", 
                min = 0, max = 1, value=.5, step =.05),
    sliderInput("pop", "Population Size:", 
                min = 10, max = 1000, value = 100, step = 1),
    sliderInput("gen", "Generations to simualte:", 
                min = 10, max = 1000, value = 100, step = 1),
    sliderInput("fit.AA", "Fitness of AA:", 
                min = 0, max = 2, value = 1, step = .05),
    sliderInput("fit.Aa", "Fitness of Aa:", 
                min = 0, max = 2, value = 1, step = .05),
    sliderInput("fit.aa", "Fitness of aa:", 
                min = 0, max = 2, value = 1, step = .05),
    sliderInput("iter", "Iterations:", 
                min = 1, max = 100, value = 10, step = 1),
    sliderInput("width", "Line width:", 
                min = .2, max = 6, value = 2, step = .1),
    selectInput("var.plot", "Plot:",
                list("AA" = 1, 
                     "Aa" = 2, 
                     "aa" = 3,
                     "A"  = 4,
                     "a"  = 5)),
    checkboxInput(inputId = "traj",
                  label = "Show expected outcome",
                  value = FALSE),
    actionButton("seed.val", 'Simulate')
  ),  
  mainPanel(
       h4(textOutput("caption1")),
       h4(textOutput("caption2")),
       plotOutput("genePlot", width='100%', height='600px')
  )
))