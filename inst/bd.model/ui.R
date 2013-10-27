# Define UI for bd trees
shinyUI(pageWithSidebar(
  headerPanel("Randomly Generated BD trees"),    # Application title
  sidebarPanel(
    sliderInput("num.taxa", "Number of Taxa:", 
                min = 4, max = 100, value=.5, step = 1),
    sliderInput("birth", "Number of Taxa:", 
                min = .1, max = 1, value=.5, step = .1),
    sliderInput("death", "Number of Taxa:", 
                min = 0, max = 1, value=0, step = .1),
    checkboxInput(inputId = "extinct",
                  label = "Prune Extinct Taxa",
                  value = FALSE),
    actionButton("seed.val", 'Simulate')
  ),  
  mainPanel(
    plotOutput("treePlot", width='100%', height='600px')
  )
))