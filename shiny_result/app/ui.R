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
      choices = c(run.list),
      selected = run.default
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
    ),
    HTML(
      '<div id="footer">
      <p>For more information:
      <ul>
      <li><a href="https://github.com/gladstone-institutes/Interactive-Enrichment-Analysis/blob/main/README.md"
      target="_blank">User guide</a></li>
      <li><a href="https://github.com/gladstone-institutes/Interactive-Enrichment-Analysis"
      target="_blank">Code repo</a></li>
      <li><a href="https://github.com/gladstone-institutes/Interactive-Enrichment-Analysis/blob/main/README.md#how-to-cite"
      target="_blank">How to cite</a></li>
      </ul></p>
      <p>Built by <a href="https://gladstone.org/science/bioinformatics-core"
      target="_blank">Gladstone Bioinformatics Core</a></p>
      <p>Powered by <a href="http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html"
      target="_blank">clusterProfiler</a>, 
      <a href="https://bioconductor.org/packages/release/bioc/html/enrichplot.html"
      target="_blank">enrichPlot</a>,
      <a href="https://CRAN.R-project.org/package=shiny"
      target="_blank">R Shiny</a>,
      <a href="https://new.wikipathways.org"
      target="_blank">WikiPathways</a>,
      <a href="https://pfocrwikipathways.org"
      target="_blank">Pathway Figures OCR</a>,
      <a href="https://geneontology.org/"
      target="_blank">Gene Ontology</a>
      </p></div'
    )
  ),
  
  #body content
  dashboardBody(
    #js head injection
    tags$link(
      rel="stylesheet",
      type = "text/css",
      href="https://cdn.drugst.one/libs/drugstone-buttons/0.0.1/drugstone-buttons.min.css"
    ),
    # sidebar footer
    tags$head(
      tags$style(HTML("
      #sidebarItemExpanded > :last-child {
        position: absolute;
        bottom: 0;
        width: 100%;
        background-color: #555555;
      }

    "))),
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
          HTML('<a style="float:right;font-size:1.8em;margin-top:-20px;color:darkcyan;"
        href="https://github.com/gladstone-institutes/Interactive-Enrichment-Analysis#data-tab"
        title="Learn more about this tab" target="_blank">&emsp;&#9432;</a>'),
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
              fluidRow(
                column(width=5,
                  selectInput(
                    "plot0",
                    "Choose a plot type for data:",
                    choices = c("Volcano plot",
                                "Heatmap",
                                "Bar plot")
                  )
                ),
                column(width = 3, style = "margin-top: 25px;",
                       downloadButton("download.plot.data", "Download")
                )),
              plotOutput("plot.data")
            )
          )
        ),
        #RESULT TAB
        tabPanel(
          "Results",
          HTML('<a style="float:right;font-size:1.8em;margin-top:-20px;color:darkcyan;"
        href="https://github.com/gladstone-institutes/Interactive-Enrichment-Analysis#results-tab"
        title="Learn more about this tab" target="_blank">&emsp;&#9432;</a>'),
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
              fluidRow(
                column(width=5,
                       selectInput(
                         "plot1",
                         "Choose a plot type for top results:",
                         choices = c("Dot plot (gene ratio)",
                                     "Dot plot (count)",
                                     "Emap plot",
                                     "Concept network",
                                     "Heatmap (GSEA)",
                                     "Heatmap (ORA)",
                                     "Upset (ORA)")
                       )
                ),
                column(width = 3, style = "margin-top: 25px;",
                       downloadButton("download.plot1.result", "Download")
                ))
            ),
            column(
              width = 5,
              fluidRow(
                column(width=5,
                       selectInput(
                         "plot2",
                         "Choose a plot type for a selected result:",
                         choices =c("GSEA score","STRING network")
                       )
                ),
                column(width = 3, style = "margin-top: 25px;",
                       downloadButton("download.plot2.result", "Download")
                ))
            )
          ),
          fluidRow(
            #result plots
            column(
              width = 7,
              plotOutput("plot1.result"),
              wellPanel("Plot options",
                        numericInput("showCategory","Top n results",
                                     20,min=2, max=40,step=1, width = 60),
                        numericInput("geneNumber","Top n genes",
                                     30,min=2, max=60,step=1, width = 60)
              )
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