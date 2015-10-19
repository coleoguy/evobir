library(shiny)
# Define UI for application
shinyUI(fluidPage(
  # Application title
  titlePanel("Brownian Motion"),
  # Sidebar for user input
  sidebarLayout(
    sidebarPanel(
      # select how many generations to run
      sliderInput("gens",
                  "Generations:",
                  min = 10,
                  max = 500,
                  value = 500),
      # select the number of replicates
      sliderInput("reps",
                  "Replicates:",
                  min = 10,
                  max = 500,
                  value = 100),
      # select the rate of evolution
      sliderInput("rate",
                  "Rate of evolution:",
                  min = 0,
                  max = 10,
                  value = 1),
      # select the part to plot a histogram of
      sliderInput("mon",
                  "Generation to plot:",
                  min = 10,
                  max = 500,
                  value = 100),
      # refresh via a seed change
      actionButton("seed.val", 'Refresh')
    ),
    # Show a plot of the simulation result and make it a bit larger with height argument
    mainPanel(
      plotOutput("distPlot", height="700px")
    )
  )
))
