# Define UI for bd trees
shinyUI(pageWithSidebar(
  headerPanel("Pairwise measures"),    # Application title
  sidebarPanel(
      sliderInput("radius", "radius:", 
                  min = 0, max = 100, value=40, step = .1),
      sliderInput("individuals", "individuals:", 
                  min = 2, max = 100, value=5, step = 1),
      sliderInput("iterations", "iterations:", 
                  min = 2, max = 1000, value=1000, step = 1),
      sliderInput("observed", "observed:", 
                  min = 0, max = 100, value=20, step = .1),
      submitButton("Submit"),
  tags$div(class="header", checked=NA,
           tags$p("To read a blog post about this."),
           tags$a(href="http://coleoguy.blogspot.com/2015/04/beetle-behavior-do-they-avoid-aggregate.html", "Click Here!")
)),
  mainPanel(
    plotOutput("densPlot", width='100%', height='500px')
  )
))