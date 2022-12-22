library(shinydashboard)
library(shiny)
library(shinyjs)
library(shinyBS)
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
    useShinyjs(),
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
      .wellbtn {
        background:#FFFFFF;
      }
      .plainUl {
        padding:0px;margin:0px
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
                                "Bar plot")
                  )
                )),
              plotOutput("plot.data"),
              wellPanel("Plot options", 
                        flowLayout(cellArgs = list(
                          style = "margin: 0px;width:auto"),
                          numericInput("topLab","Top n genes",
                                       10, min=0, max=20, step=1, width = 80),
                          p("or"),
                          selectizeInput(
                            inputId = "selectLab",
                            label = "Enter genes by name",
                            choices = NULL,
                            multiple = TRUE,
                            width = "100%",
                            options = list(
                              'plugins' = list('remove_button'),
                              'create' = TRUE
                            )
                          ),
                          br(),
                          numericInput("category_label_0","Label font size",
                                       12,min=9, max=24,step=1, width = 90),
                          selectInput("legend_pos","Legend position",
                                      choices = c("top", "right", "bottom", "left"),
                                      selected = "right", width = 110)
                        )
              ),
              wellPanel("PDF options",
                        fluidRow(
                          column(width = 3, style = "margin-top: 5px;",
                                 downloadButton("download.plot.data", "Download", class = "wellbtn")
                          ),
                          column(width = 3,style = "margin-top: -20px;padding: 0 5px 0 5px;",
                                 numericInput("plot.width", "width (px)", 1200, 
                                              min = 200, max=1200, step = 100)
                          ),
                          column(width = 3,style = "margin-top: -20px;padding: 0 5px 0 5px;",
                                 numericInput("plot.height", "height (px)", 1200, 
                                              min = 200, max=1200, step = 100)
                          )
                        ))
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
          br(),
          fluidRow(
            #selections and download buttons for plots
            column(width=5,
                   selectInput(
                     "plot1",
                     "Choose a plot type for top results:",
                     choices = c("Dot plot",
                                 "Bar plot",
                                 "Heatmap",
                                 "Emap plot",
                                 "Concept network")
                   )
            ),
            column(width=5, offset=2,
                   selectInput(
                     "plot2",
                     "Choose a plot type for a selected result:",
                     choices =c("GSEA score")
                   )
            )
          ),
          fluidRow(
            #result plots
            column(
              width = 7,
              plotOutput("plot1.result"),
              wellPanel("Plot options", 
                        flowLayout(cellArgs = list(
                          style = "margin: 0px;width:auto"),
                        numericInput("showCategory","Top n results",
                                     20, min=2, max=40, step=1, width = 80),
                        numericInput("showGene","Top n genes",
                                     30, min=2, max=60, step=1, width = 80),
                        selectInput("dotplot_x","X-axis by",
                                    choices = c("Count","GeneRatio"), 
                                    selected = "Count", width=110),
                        shinyBS::bsTooltip("dotplot_x", 
                                           paste('<ul class="plainUl">',
                                                 "<li>Count = number of dataset genes overlapping a given term</li>",
                                                 "<li>GeneRatio = dataset genes overlapping a given term vs all possible terms</li>",
                                                 "</ul>"), 
                                           placement = "right", trigger = "hover"),
                        selectInput("dotplot_size","Size by",
                                    choices = c("Count","GeneRatio","Percentage"), 
                                    selected = "Percentage", width=110),
                        shinyBS::bsTooltip("dotplot_size", 
                                           paste('<ul class="plainUl">',
                                                 "<li>Count = number of dataset genes overlapping a given term</li>",
                                                 "<li>GeneRatio = dataset genes overlapping a given term vs all possible terms</li>",
                                                 "<li>Percentage = dataset genes overlapping a given term vs all genes in term</li>",
                                                 "</ul>"), 
                                           placement = "right", trigger = "hover"),
                        selectInput("bar_by","X-axis by",
                                    choices = c("Fold change","Count","Percentage"), 
                                    selected = "Fold change", width=120),
                        shinyBS::bsTooltip("bar_by", 
                                           paste('<ul class="plainUl">',
                                                 "<li>Fold change = average fold change of overlapping genes</li>",
                                                 "<li>Count = number of dataset genes overlapping a given term</li>",
                                                 "<li>Percentage = dataset genes overlapping a given term vs all genes in term</li>",
                                                 "</ul>"), 
                                           placement = "right", trigger = "hover"),
                        numericInput("label_format","Wrap labels",
                                     50, min=10, max=80, step=1, width = 80),
                        numericInput("category_label","Label font size",
                                     12,min=8, max=24,step=1, width = 90),
                        numericInput("gene_label","Gene font size",
                                     12,min=8, max=24,step=1, width = 90),
                        numericInput("category_node","Result node scale",
                                     1,min=0.5, max=4,step=0.5, width = 110),
                        numericInput("gene_node","Gene node scale",
                                     1,min=1, max=12,step=2, width = 110),
                        selectInput("color","Color by",
                                    choices = c("p.adjust","pvalue","qvalue"),
                                    selected = "p.adjust", width = 90),
                        selectInput("layout","Layout",
                                    choices = c("nicely","kk","sugiyama","fr", "gem","lgl","mds"),
                                    selected = "nicely", width = 105),
                        selectInput("cnet_layout","Layout",
                                    choices = c('star', 'circle', 'gem', 'dh', 'graphopt',  
                                                'mds', 'fr', 'kk', 'drl', 'lgl'),
                                    selected = "kk", width = 105),
                        checkboxInput("cnet_circular","Circular",
                                    value = FALSE, width = 80),
                        checkboxInput("cnet_colorEdge","colorEdge",
                                    value =  TRUE, width = 80)
                        )
              ),
              wellPanel("PDF options",
                        fluidRow(
                          column(width = 3, style = "margin-top: 5px;",
                                 downloadButton("download.plot1.result", "Download", class = "wellbtn")
                          ),
                          column(width = 2,style = "margin-top: -20px;padding: 0 5px 0 5px;",
                                 numericInput("plot1.width", "width (px)", 1200, 
                                              min = 200, max=1200, step = 100)
                          ),
                          column(width = 2,style = "margin-top: -20px;padding: 0 5px 0 5px;",
                                 numericInput("plot1.height", "height (px)", 1200, 
                                              min = 200, max=1200, step = 100)
                          )
                        ))
            ),
            column(
              width = 5,
              htmlOutput("html.result"),
              plotOutput("plot2.result"),
              conditionalPanel(
                condition = "input.plot2 == 'GSEA score'",
                wellPanel("PDF options",
                          fluidRow(
                            column(width = 4, style = "margin-top: 5px;",
                                   downloadButton("download.plot2.result", "Download", class = "wellbtn")
                            ),
                            column(width = 3,style = "margin-top: -20px;padding: 0 5px 0 5px;",
                                   numericInput("plot2.width", "width (px)", 400, 
                                                min = 100, max=600, step = 100)
                            ),
                            column(width = 3,style = "margin-top: -20px;padding: 0 5px 0 5px;",
                                   numericInput("plot2.height", "height (px)", 300, 
                                                min = 100, max=600, step = 100)
                            )
                          ))
              )
            )
          )
        ) #end tabPanel
      )
    )
  )
)