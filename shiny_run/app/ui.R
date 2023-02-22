library(shinydashboard)
library(shiny)
library(shinyjs)
library(shinyBS) #see bsTooltip in server.R
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
    uiOutput('run.button'),
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
    useShinyjs(),
    extendShinyjs(text = jscode, functions = c("togglebox", #see global.R for jscode 
                                               "hidebox",
                                               "showbox",
                                               "startnew")), 
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
      #debug
      textOutput('debug.text')
    ),
    fluidRow(
      box(id = "db",
        title = "Databases", status = "warning", solidHeader = TRUE, collapsible = T,
        HTML('<a style="float:right;font-size:1.8em;margin-top:-10px;color:darkorange;"
        href="https://github.com/gladstone-institutes/Interactive-Enrichment-Analysis#databases"
        title="Learn more about this panel" target="_blank">&emsp;&#9432;</a>'),
        selectInput(
          "rdata",
          "Choose a database collection",
          choices = c("",rdata.list),
          selected = NULL,
          multiple = F
        ),
        uiOutput("db.build"),
        textOutput("db.list.title"),
        tableOutput("db.list"),
        uiOutput("gmt.upload"),
        htmlOutput("db.ready")
      ),
      box(id = "ds",
        title = "Datasets", status = "info", solidHeader = TRUE, collapsible = T,
        HTML('<a style="float:right;font-size:1.8em;margin-top:-10px;color:darkskyblue;"
        href="https://github.com/gladstone-institutes/Interactive-Enrichment-Analysis#datasets"
        title="Learn more about this panel" target="_blank">&emsp;&#9432;</a>'),
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
        uiOutput("plan.params"),
        htmlOutput("ds.ready")
      ),
      uiOutput("progress")
    )
  )
)