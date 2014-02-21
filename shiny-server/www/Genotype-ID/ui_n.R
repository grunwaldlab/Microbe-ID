library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Input Table"),
  sidebarPanel(    textInput("mst1","","100/100"),
                   textInput("mst2","","100/100"),
                   textInput("mst3","","100/100"),
                   textInput("mst4","","100/100"),
                   textInput("mst5","","100/100"),
                   textInput("mst6","","100/100"),
                   textInput("mst7","","100/100"),
                   textInput("mst8","","100/100"),
                   textInput("mst9","","100/100")
                   ,
                   submitButton("Calculate Tree"),
                   downloadButton('downloadData', 'Download NEWICK Tree')
  ),
  mainPanel(
    plotOutput("distPlotTree")
    )
))
