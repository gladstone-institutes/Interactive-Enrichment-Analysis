#load required packages
require(EnhancedVolcano)
require(writexl)
require(clusterProfiler)
require(ggplot2)
require(enrichplot)

shinyServer(function(input, output, session) {
  
  rv <- reactiveValues()
  
  ########### 
  # sidebar #
  ###########
  
  #############
  # databases #
  #############
  
  getDatabaseSet <- function(){
    db.list = NULL
    output$db.list.title <- NULL
    rv$db.status <- FALSE
    if (!input$rdata == ""){
      if(input$rdata == "BUILD NEW DATABASE"){
        #TODO list GMT files... 
      } else {
        rdata.fn <- file.path("../databases",input$rdata)
        db.list <- load(rdata.fn)
        for (db in db.list){
          db.colnames <- tolower(names(eval(parse(text=db))))
          for (cn in c("name","term","gene")){
            if(!cn %in% db.colnames){
              stop(sprintf('%s is missing a %s column',db, cn))
            }
          }
        }
        output$db.list.title <- renderText("These are the databases we will use in this run:")
        rv$db.status <- TRUE
      }
    }
    return(db.list)
  }
  
  output$db.list <- renderTable({
    getDatabaseSet()
  },
  colnames = F,
  hover = T)
  
  
  ############
  # datasets #
  ############
  
  getDatasets <- function(){
    output$sample.ds.title <- NULL
    output$gene.check <- NULL
    rv$ds.status <- FALSE
    for (input.ds in input$datasets){
      if(input.ds =="ADD NEW DATASETS"){
        updateSelectInput(session, "datasets", choices = c(""), selected = c(""))
        uploadDataset()
        return()
      } else if (input.ds =="SELECT ALL"){
        updateSelectInput(session, "datasets", selected = list.files("../datasets", "csv"))
        return()
      }
    }
    if (length(input$datasets) > 0){
      # output$debug.text<-renderText(paste(input$datasets,collapse = ","))
      ds.fn <- file.path("../datasets",input$datasets[1])
      data.df <- read.table(ds.fn, sep = ",", header = T, stringsAsFactors = F)
      checkTableData(data.df)
      return(data.df)
    }
  }
  
  #check table data
  checkTableData <- function(data.df) {
    output$sample.ds.title <- renderText("Let's examine as sample of the first one...")
    ds.names <- tolower(names(data.df))
    #GENE CHECK
    if(!'gene' %in% ds.names){
      stop('Please reformat your CSV to have a "gene" column with gene names.')
    } else {
      output$gene.check <- renderUI({
        selectInput(
          "from.type",
          "Choose type of gene identifier",
          choices = supported.idTypes,
          multiple = F
        )
      })
      rv$fromType <- input$from.type
    }
    #TODO
    rv$ds.status <- TRUE
  }
  
  #render table data
  output$table.sample.ds <- DT::renderDataTable(server = TRUE,{
    DT::datatable({ head(getDatasets()) },
                  rownames=FALSE,
                  options = list(
                    dom = 'Brt'
                  )
    )
  })
  
  #upload dataset csv
  uploadDataset <- function(){
    output$sample.ds.title <- NULL
    output$ds.upload <- renderUI({
      column( width = 12,
      fileInput("ds.file", "Upload CSV Files",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      actionButton("cancel.upload", "Cancel")
      )
    })
  }
  
  #when upload is cancels
  observeEvent(input$cancel.upload, {
    updateSelectInput(session, "datasets", selected = c(input$ds.file$name),
                      choices = c(list.files("../datasets", "csv"),
                                  "SELECT ALL", "ADD NEW DATASETS"))
    output$ds.upload <- NULL
    getDatasets()
  })
  
  #when upload is performed
  observeEvent(input$ds.file, { 
    # output$debug.text<-renderText(paste(input$ds.file$name,collapse = ","))
    uploaded.files = file.rename(input$ds.file$datapath, paste0(input$ds.file$name))
    uploaded.files = paste0(input$ds.file$name)
    file.copy(uploaded.files, "../datasets/.")
    file.remove(uploaded.files)
    updateSelectInput(session, "datasets", selected = c(input$ds.file$name),
                      choices = c(list.files("../datasets", "csv"),
                                  "SELECT ALL", "ADD NEW DATASETS"))
    output$ds.upload <- NULL
    getDatasets()
  })
  
  #set run status
  observe({
    if(rv$ds.status){
      output$ds.ready <-  renderText("&emsp;&#10004; Ready!")
      output$ds.status <- renderText("&emsp;&#10004; Ready!")
    } else {
      output$ds.ready <- NULL
      output$ds.status <- NULL
    }
    if(rv$db.status){
      output$db.ready <-  renderText("&emsp;&#10004; Ready!")
      output$db.status <- renderText("&emsp;&#10004; Ready!")
    }else {
      output$db.ready <- NULL
      output$db.status <- NULL
    }
    if(rv$db.status & rv$ds.status){
      output$run.status <- renderText("<h3>&nbsp;&nbsp;3. Let's go!</h3>")
      output$run.button <- renderUI({
        actionButton("run", "Run analysis")
      })
    } else {
      output$run.status <- NULL
      output$run.button <- NULL
    }
  })
  
  # 
  # #differential analysis tab - select the cluster number
  # output$DEclusterNumUI <- renderUI({
  #   selectInput(
  #     "clusterNum",
  #     "Choose a cluster to compare with all other clusters for differentially expressed genes",
  #     choices = sort(unique(get(input$dataset)@meta.data$seurat_clusters)),
  #     selected = sort(unique(get(input$dataset)@meta.data$seurat_clusters))[1]
  #   )
  # })
  
  
})
