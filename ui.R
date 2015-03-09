library(shiny)
library(ggplot2)
library(lattice)
library(RColorBrewer)
library(quantmod)
library(grid)
library(flux)
library(shinysky)
library(shinyBS)
options(shiny.maxRequestSize=30*1024^2)
options(java.parameters = "-Xmx8000m")

shinyUI(pageWithSidebar(
  headerPanel("Ribosome Profiling Analysis"),
  sidebarPanel(fileInput('file1',label=h4("Upload ribosome profiling data file:"),accept=c('csv','comma-separated-values'),multiple = TRUE),
               radioButtons("mode","",list("Quick Plot"="quick","Main Mode"="main")),
               conditionalPanel(condition="input.mode == 'main'",
                                list(
                                  checkboxInput("check",h5("Name the sample"),FALSE),
                                  conditionalPanel(condition="input.check == TRUE",
                                                   uiOutput("insertname")),
                                  checkboxInput("check2",h5("Group the sample"),FALSE),
                                  conditionalPanel(condition="input.check2 == TRUE",
                                                   uiOutput("groupingnumber"),
                                                   uiOutput("groupno"))
                                  )
                                ),
               radioButtons("colour1",h5("Total graph color collocation:"),list("Classical"="classical","Similar colors"="similar","Same color"="same")),
               conditionalPanel(condition="input.colour1 == 'similar' | input.colour1 == 'same'",
                                uiOutput("totalcolour")),
               sliderInput("lineweight",h5("Line weight:"),min=0,max=1.8,value=0.9,step=0.05),               
               busyIndicator(text= HTML(paste("Calculation in progress..","Please click again when stuck.",sep="<br/>")), wait = 0),
               actionButton("goButton", "Analyze",styleclass="primary"),
               helpText("."),
               h4("Note:"),
               helpText("This app is used to analyze ribosome profiling data files (standard csv format) created by gradiant fractionation instrument. Multiple files or files with multigroups analysis are also available. In 'Quick Plot' mode, you can just upload data and click 'Analyze' button to get a glance. Using 'Main Mode', graphs could be easily beautified by setting group names and adjusting color styles or line weight. After inputing data, real name of sample can be typed in by 'Name the sample' option. Clicking 'Group the sample', typing the exact number and name of group, then you can select group corresponding to each sample. Considering about aesthetics, 2 groupings and 6 groups are the top numbers we can support. Seperated AUC graphs of polysome or monosome are only visible in no-grouping mode. App would have a delay during graphs generation or renewal, please be patient."),
               helpText("Attention: The first absorbance value of each sample should below 0.03, otherwise some unknown errors will occur. Large data analysis needs more time and may lead to program breakdown. Weird curve trend may lead to incorrect normalization or calculation. In case of app crash or errors, we just need to refresh it.")
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel(width=18,
               h4("Overview"),
               uiOutput("singlecolor"),
               plotOutput("graph1"),
               uiOutput("NormalizeDIY"),
               uiOutput("DIYModal"),
               plotOutput("graph2")),
      tabPanel(width=18,
               h4("Seperated Curve"),
               uiOutput("singlesetting"),
               div(class="row",div(class="span1",""),div(class="span6",shinyalert("shinyalert1"))),
               uiOutput("lattice"),               
               uiOutput("plotchoose"),
               plotOutput("singleplot")),
      tabPanel(h4("AUC Ratio"),
               uiOutput("Showauc"),
               plotOutput("graph3"),
               plotOutput("graph4"),
               plotOutput("graph6"),
               plotOutput("graph5")),
      tabPanel(h4("Raw Data"),
               h5(dataTableOutput("table"))),
      tabPanel(h4("Information"),
               tableOutput("information"))
    ))
  )
)
