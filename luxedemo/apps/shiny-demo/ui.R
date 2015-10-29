# front end of LUXE demo app

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("The LUXE Demo App using Shiny"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      helpText("This demo app creates an interactive version of Figure 1 in Huang et al. (2013)"),
      helpText("Parameters for column 1"),
      uiOutput("column1"),
      helpText("Parameters for column 2"),
      uiOutput("column2"),
      helpText("Parameters for column 3"),
      uiOutput("column3")      
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))
