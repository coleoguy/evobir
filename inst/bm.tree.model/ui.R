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
                  value = 200),
      # select the number of replicates
      sliderInput("birth",
                  "Birth Rate:",
                  min = .01,
                  max = .02,
                  value = .015),
      # select the part to plot a histogram of
      sliderInput("mon",
                  "Generation to plot:",
                  min = 10,
                  max = 500,
                  value = 200),
      # refresh via a seed change
      actionButton("seed.val", 'Refresh')
    ),
    # Show a plot of the simulation result and make it a bit larger with height argument
    mainPanel(
      plotOutput("distPlot", height="600px")
    )
  )
))
