library(shiny)
shinyUI(
  # we will use h2 for all major divisions
  fluidPage(
    titlePanel("Visualizing Trees and Discrete Data"),
    h2("1. Load data"),
    fluidRow(
      column(4,
             "Select tree",
             fileInput("tree", label = ""),
             h6("tree should be in newick format and the tip names must match the data file")),
      
      column(4,
             "Select data",
             fileInput("traits", label = ""),
             checkboxInput("headers",
                           "Headers present", value = TRUE),
             h6("data should be in a CSV file with two columns first should be species names second should be the trait of interest")),
      column(4,
             selectInput("plotType", "",
               c("tip states" = "tips",
                 "ML ancestral states" = "mlasr",
                 "Stochastic mapping" = "stmap")))
    ),
    # Create a new row for the tree.
    h2("2. Visualize data"),
    "Select the information to include with the phylogeny:",
    fluidRow(
      column(3,
        selectInput("type", h6("Plot type"), 
                    list("phylogram" = "phylogram",
                         "fan" = "fan")),
        sliderInput("cex.font", label = h6("Font size"), min = .1,
                    max = 5, step = .1, value = 1),
      
      
      # each call of conditional panel is for a different type of
      # plot.
      
      ## Stochastic Mapping
      conditionalPanel(
        condition = "input.plotType == 'stmap'",
        selectInput("modelSIMMAP", h6("Trait model"), 
                    list("different rates" = "ARD",
                         "symmetrical" = "SYM"))),
      
      
      ## ML ASR Mapping
      conditionalPanel(
        sliderInput("cex.tipML", label = h6("Symbol size"), min = .1,
                    max = 5, step = .1, value = 1),
        condition = "input.plotType == 'mlasr'",
        selectInput("modelACE", h6("Trait model"), 
                    list("different rates" = "ARD",
                         "symmetrical" = "SYM"))),
      
      ## Just tip states
      conditionalPanel(
        condition = "input.plotType == 'tips'",
        sliderInput("pch", label = h6("Trait symbol"), min = 0,
          max = 24, step = 1, value = 16),
        sliderInput("cex.tip", label = h6("Symbol size"), min = .1,
          max = 5, step = .1, value = 1)),
      
      ## Used on all
      sliderInput("lwd", label = h6("Line width"), min = .25,
                  max = 10, step = .25, value = 3),
      textInput("col1", h6("State 1 color"), value = "#ffcc33"),
      textInput("col2", h6("State 2 color"), value = "#7a0019")      
      
      
      
      
    ),
    
    
    
    column(9, plotOutput("plot.smap", height = "600px"))
    )
  )
)
