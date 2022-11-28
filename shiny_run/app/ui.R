#load required packages
library(shinydashboard)
library(shiny)
library(shinyjs)
library(shinyBS)
library(DT)

dashboardPage(
  
  #header
  dashboardHeader(
    #application title
    title = "Enrichment Analysis Run"
  ),
  
  #sidebar content
  dashboardSidebar(
    h3(HTML('&nbsp;'),'1. Databases'),
    htmlOutput('db.status'),
    h3(HTML('&nbsp;'),'2. Datasets'),
    htmlOutput('ds.status'),
    htmlOutput('run.header'),
    uiOutput('run.button')
  ),
  
  #body content
  dashboardBody(
    useShinyjs(),
    extendShinyjs(text = jscode, functions = c("togglebox", #see global.R for jscode 
                                               "hidebox",
                                               "showbox",
                                               "startnew")), 
    fluidRow(
      #debug
      textOutput('debug.text')
    ),
    fluidRow(
      box(id = "db",
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
      box(id = "ds",
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
        htmlOutput("ds.msg"),
        uiOutput("ds.upload"),
        br(),
        uiOutput("run.params"),
        uiOutput("ora.params"),
        htmlOutput("analysis.plan"),
        htmlOutput("ds.ready")
      ),
      uiOutput("progress")
    )
  )
)