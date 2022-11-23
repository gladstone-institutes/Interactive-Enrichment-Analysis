#load required packages
library(shinydashboard)
library(shiny)
library(DT)

dashboardPage(
  
  #header
  dashboardHeader(
    #application title
    title = "Enrichment Analysis Run"
  ),
  
  #sidebar content
  dashboardSidebar(
    h3(' 1. Databases'),
    htmlOutput('db.status'),
    h3(' 2. Datasets'),
    htmlOutput('ds.status'),
    htmlOutput('run.status'),
    uiOutput('run.button')
  ),
  
  #body content
  dashboardBody(
    fluidRow(
      #debug
      textOutput('debug.text')
    ),
    fluidRow(
      box(
        title = "Databases", status = "warning", solidHeader = TRUE, collapsible = T,
        selectInput(
          "rdata",
          "Choose a database collection",
          choices = c("",rdata.list),
          selected = NULL,
          multiple = F
        ),
        textOutput("db.list.title"),
        tableOutput("db.list"),
        htmlOutput("db.ready")
      ),
      box(
        title = "Datasets", status = "info", solidHeader = TRUE, collapsible = T,
        selectInput(
          "datasets",
          "Choose one or more datasets",
          choices = c("",ds.list),
          selected = NULL,
          multiple = T
        ),
        textOutput("sample.ds.title"),
        div(DT::dataTableOutput("table.sample.ds"), style = "font-size:80%; line-height:30%"),
        uiOutput("ds.upload"),
        br(),
        uiOutput("gene.check"),
        htmlOutput("ds.ready")
      )
    )
  )
)