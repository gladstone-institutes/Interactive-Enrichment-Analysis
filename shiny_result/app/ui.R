#load required packages
library(shinydashboard)
library(shiny)
library(DT)

dashboardPage(
  
  #header
  dashboardHeader(
    #application title
    title = "Enrichment Analysis Explorer"
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
        includeHTML("../summary.txt")
      )
    ),
    fluidRow(
      #debug
      textOutput('debug.text')
    ),
    fluidRow(
      tabBox(
        width = 12,
        #DATA TAB
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
              div(DT::dataTableOutput("table.data"), style = "font-size:80%; line-height:30%")
              # DT::dataTableOutput("table.data")
            ),
            #data plot
            column(
              width = 6,
              downloadButton("download.plot.data", "Download"),
              plotOutput("plot.data")
            )
          )
        ),
        #RESULT TAB
        tabPanel(
          "Results",
          selectInput(
            "database",
            "Choose a database",
            choices = c(db.list)
          ) ,
          #results table
          strong("Select a table row to update plot data"),
          fluidRow(
            column(
              width = 12,
              div(DT::dataTableOutput("table.result"), style = "font-size:80%; line-height:30%")
            )
          ) ,
          #result plots
          fluidRow(
            #selections and download buttons for plots
            column(
              width = 7,
              selectInput(
                "plot1",
                "Choose a plot type for top results:",
                choices = c("Dot plot (gene ratio)",
                            "Dot plot (count)",
                            "Emap plot",
                            "Concept network",
                            "Concept network (circular)",
                            "Volcano plot (GSEA)",
                            "Volcano plot (ORA)",
                            "Heatmap (GSEA)",
                            "Heatmap (ORA)",
                            "Upset (ORA)")
              ),
              downloadButton("download.plot1.result", "Download")
            ),
            column(
              width = 5,
              selectInput(
                "plot2",
                "Choose a plot type for a selected result:",
                choices =c("GSEA score","STRING network")
              )
            ),
            downloadButton("download.plot2.result", "Download")
          ),
          fluidRow(
            #result plots
            column(
              width = 7,
              plotOutput("plot1.result")
            ),
            column(
              width = 5,
              htmlOutput("html.result"),
              plotOutput("plot2.result")
            )
          )
        ) #end tabPanel
      )
    )
  )
)