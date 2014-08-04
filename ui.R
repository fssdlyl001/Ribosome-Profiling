library(shiny)
library(ggplot2)
library(quantmod)
library(grid)
options(shiny.maxRequestSize=30*1024^2)
options(java.parameters = "-Xmx8000m")

shinyUI(pageWithSidebar(
  headerPanel("Ribosome Profiling Analysis"),
  sidebarPanel(fileInput('file1',label=h4("Upload ribosome profiling data file:"),accept=c('csv','comma-separated-values'),multiple = TRUE),
               radioButtons("colour1",h5("Total graph color collocation:"),list("Classical"="classical","Similar Colors"="similar")),
               conditionalPanel(condition="input.colour1 == 'similar'",
                                selectInput('similar1','Select hue:',choices=c("Cold" = "cold", "Warm" = "warm", "Grey" = "grey"))),
               selectInput("colour2", h5("Single graphs color:"),choices=c("Blue" = "blue", "Red" = "red", "Orange" = "orange", "Green"="green", "Black" = "black")),
               sliderInput("lineweight",h5("Line weight:"),min=0,max=1.5,value=0.75,step=0.05),
               actionButton("goButton", "Analyze"),
               helpText("."),
               checkboxInput("check",h5("Name the group"),FALSE),
               conditionalPanel(condition="input.check == TRUE",
                                uiOutput("insertname")),
               h4("Note:"),
               helpText("This app is used to analyze ribosome profiling data files (standard csv format) created by gradiant fractionation instrument. Multiple files or files with multigroups analysis are also available. Upload your data and click 'Analyze' button, and other processes will be automatically completed by the program. In the left sidebar panel, you can type the exact sample names and easily adjust color styles or line weight. App would have a delay during graphs generation or renewal, please be patient."),
               helpText("Attention: The first absorbance value of each sample should below 0.03, otherwise some unknown errors will occur. Large data analysis needs more time and may lead to program breakdown, especially in 'Raw Data' part. In case of app crash or errors, we just need to refresh it.")
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel(width=18,
               h4("Overview"),
               plotOutput("graph1"),
               plotOutput("graph2")),
      tabPanel(width=18,
               h4("Single Graphs"),
               uiOutput("plots")),
      tabPanel(h4("Raw Data"),
               h5(tableOutput("table"))),
      tabPanel(h4("Information"),
               tableOutput("information"))
    ))
))
