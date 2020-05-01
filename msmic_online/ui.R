##ui.R
library("shiny")
library("foreign")
library("ggplot2")

shinyUI(pageWithSidebar(
    
  
  
    # Header:
    headerPanel("Mass Spectrometry Microscopy "),
    
    # Input in sidepanel:
    sidebarPanel(
        tags$style(type='text/css', ".well { max-width: 20em; }"),
        # Tags:
        tags$head(
            tags$style(type="text/css", "select[multiple] { width: 100%; height:10em}"),
            tags$style(type="text/css", "select { width: 100%}"),
            tags$style(type="text/css", "input { width: 19em; max-width:100%}")
        ),
        
        
        radioButtons("LC", "Choose localization: ", choices=c("Whole cell", "Mitochondrial", "test")),
        fileInput("datafile", "Upload data-file:"),
        htmlOutput("ArgSelect"),
        
        br(),
        p("Input file should be comma separated SAINT output file and contain at least columns \"Bait\", \"Prey\", and \"AvgSpec\"."),
        br(),
        p(strong("Remember to filter the data.")),
        downloadLink('downloadFigure', 'Download Figure as .pdf'),
        br(),
        br(),
        br(),
        br(),
        downloadLink('downloadExample1', 'Download Example input file 1'),
        br(),
        downloadLink('downloadExample2', 'Download Example input file 2'),
        br(),
        downloadLink('downloadOutput', 'Download Example output')
        
        
    ),
    
    
    mainPanel(
        plotOutput("plot"),
        tabsetPanel(
            id = 'dataset',
            tabPanel('Matrix', tableOutput('table')),
            tags$style(type="text/css",
                       ".shiny-output-error { visibility: hidden; }",
                       ".shiny-output-error:before { visibility: hidden; }"
                       
            )
        
        
    )
)))

