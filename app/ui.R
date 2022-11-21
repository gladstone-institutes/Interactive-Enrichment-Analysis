#load required packages
library(shinydashboard)
library(shiny)
library(DT)

dashboardPage(
  
  #header
  dashboardHeader(
    # Application title
    title = "Enrihcment Analysis Explorer"
  ),
  
  #sidebar content
  dashboardSidebar(
    #data set drop down menu
    selectInput(
      "run",
      "Choose a run",
      choices = c(run.list)
    ),
    selectInput(
      "dataset",
      "Choose a dataset",
      choices = c(ds.list)
    ),
    selectInput(
      "method",
      "Choose a method",
      choices = c(method.list)
    )
  ),
  
  #body content
  dashboardBody(
    fluidRow(
      #collapsible summary box with html 
      box(
        id = "sumary",
        title = "Summary",
        width = 12,
        collapsible = TRUE,
        collapsed = TRUE,
        includeHTML("summary.txt")
      )
    ),
    fluidRow(
      #debug
      textOutput('debug.text')
    ),
    fluidRow(
      #tabbed Box for plots and DE analysis
      tabBox(
        width = 12,
        #data tab
        tabPanel(
          "Data",
          #checkbox to select data view
          checkboxInput("excluded.data", 
                        "View excluded data rows:",
                        value = FALSE, 
                        width = NULL
          ),
          fluidRow(
            #data table
            column(
              width = 6,
              div(dataTableOutput("table.data"), style = "font-size:80%; line-height:30%")
              # DT::dataTableOutput("table.data")
            ),
            column(
              width = 6,
              downloadButton("download.plot.data", "Download"),
              plotOutput("plot.data")
            )
          )
        ),
        do.call(tabsetPanel, c(id='t',lapply(1:5, function(i) {
          tabPanel(
            title=paste0('tab', i),
            textOutput(paste0('a',i))
          )
        })))
        
        # lapply(1:5, function(i) {
        #   tabPanel(title = i,
        #            textOutput(paste0('a',i)))
        # })
        
        # #GO tab
        # tabPanel(
        #   "GO",
        #   h4("Enriched Gene Ontology terms: Select a row to update plots below"),
        #   fluidRow(#results table
        #     column(
        #       width = 12,
        #       DT::dataTableOutput("table.result")
        #     )
        #   ),
        #   #result plots
        #   fluidRow(
        #     #selections and download buttons for plots 
        #     column(
        #       width = 6, 
        #       selectInput(
        #         "plot1",
        #         choices = c("one","two")
        #       ),
        #       downloadButton("download.plot1.result", "Download")
        #     ),
        #     column(
        #       width = 6, 
        #       selectInput(
        #         "plot2",
        #         choices =c("one","two")
        #       ),
        #       downloadButton("download.plot2.result", "Download")
        #     )
        #   ),
        #   fluidRow(
        #     #result plots 
        #     column(
        #       width = 6, 
        #       plotOutput("plo1.result")
        #     ),
        #     column(
        #       width = 6, 
        #       plotOutput("plot2.result")
        #     )
        #   )
        # ) #end tabPanel
      )
    )
  )
)