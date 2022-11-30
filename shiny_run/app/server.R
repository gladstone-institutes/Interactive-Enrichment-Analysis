library(shinyjs)
library(rstudioapi)

shinyServer(function(input, output, session) {
  
  # output$debug.text<-renderText(rv$params$minGSSize)
  
  params = data.frame(
    db.name = "NA",
    db.list = "NA",
    org.db.name = "org.Hs.eg.db",
    fromType = "SYMBOL",
    minGSSize = minGSSize.default,
    maxGSSize = maxGSSize.default,
    ora.fc = fc.default,
    ora.pv = pv.default,
    run.ora = FALSE,
    run.gsea = FALSE
  )
  rv <- reactiveValues(params = params,
                       lib.status = FALSE,
                       logfile = "",
                       plan = "")
  
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
        db.list <- load(rdata.fn, globalenv())
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
        rv$params$db.name <- input$rdata
        rv$params$db.list <- paste(db.list, collapse = ",")
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
    output$run.params <- NULL
    output$ora.params <- NULL
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
    output$sample.ds.title <- renderText("Let's examine as sample of your dataset...")
    ds.names <- tolower(names(data.df))
    ds.names.in <- intersect(analyzed.columns, ds.names)
    ds.names.out <- setdiff(analyzed.columns, ds.names)
    #GENE COLUMN CHECK
    if(!'gene' %in% ds.names.in){
      output$ds.msg <- renderText({
        paste("<p><font color=\"#FF0000\"><b>",
        "Error: dataset missing <i>gene</i> column.</b>",
        "<br />Please reformat your CSV to have a column named <i>gene</i> 
        with symbols or identifiers.",
        "</font>")
      })
      rv$params$run.ora <- FALSE
      rv$params$run.gsea <- FALSE
    } else {
      output.ds.msg.optional <- ""
      if(length(ds.names.out) > 0)
        output.ds.msg.optional <- paste(
                               "<br />Optional columns not found:",
                               "<i>",paste(ds.names.out, collapse =", "),"</i>")
      output.ds.msg <-paste("<p><b>",
                            "Required columns found:",
                            "<i>",paste(ds.names.in, collapse =", "),"</i>",
                            "</b>",output.ds.msg.optional,"</p>")
      output$ds.msg <- renderText({ output.ds.msg })
      #SET METHODS and PLAN
      rv$params$run.ora <- TRUE
      rv$plan <- "Analysis plan: <b>run ORA only</b> (GSEA requires <i>rank</i> or <i>p.value</i> column)"
      if('p.value' %in% ds.names.in | 'rank' %in% ds.names.in ){
        rv$params$run.gsea <- TRUE
        rv$plan <- "Analysis plan: <b>run ORA and GSEA</b>"
      } else {
        rv$params$run.gsea <- FALSE
      }
      #SET PARAMS
      #choose id type (fromType for ID mapping)
      output$run.params <- renderUI({
        tagList(
          selectInput(
            "fromtype",
            "Choose type of gene identifier",
            choices = supported.idTypes,
            multiple = F,
            selected = rv$params$fromType
          ),
          shinyBS::bsTooltip("fromtype", 
                    "Which type of identifier do you have in the <i>gene</i> column above?", 
                    placement = "bottom", trigger = "hover"),
          #choose org name (org.db.name for ID mapping)
          selectInput(
            "orgname",
            "Choose organism",
            choices = supported.orgs,
            multiple = F,
            selected = rv$params$org.db.name
          ),
          shinyBS::bsTooltip("orgname", 
                    "Which organism is represented by the gene identifiers above?", 
                    placement = "bottom", trigger = "hover"),
          #set min/maxGSSize
          fluidRow(
            column( 
              width = 6,
              numericInput(
                "minGSSize",
                "Set minimum gene set size",
                value = rv$params$minGSSize,
                min = 1, max = 50,
                step = 1
              ),
              shinyBS::bsTooltip("minGSSize", 
                        "The minimum gene set size to be included in the analysis, i.e., minGSSize", 
                        placement = "bottom", trigger = "hover")
            ),
            column( 
              width = 6,
              numericInput(
                "maxGSSize",
                "Set maximum gene set size",
                value = rv$params$maxGSSize,
                min = 100, max = 1000,
                step = 100
              ),
              shinyBS::bsTooltip("maxGSSize", 
                        "The maximum gene set size to be included in the analysis, i.e., maxGSSize", 
                        placement = "bottom", trigger = "hover")
            )
          )
        )
      })
      # ORA ONLY: fold.change and p.value
      #set min/maxGSSize
      if ('fold.change' %in% ds.names & 'p.value' %in% ds.names){
        output$ora.params <- renderUI({
          tagList(
            fluidRow(
              column( 
                width = 6,
                numericInput(
                  "foldchange",
                  "Set fold change cutoff",
                  value = rv$params$ora.fc,
                  min = 0, max = 5,
                  step = 0.1
                ),
                shinyBS::bsTooltip("foldchange", 
                          "Absolute value threshold for defining subset of genes for ORA analysis", 
                          placement = "bottom", trigger = "hover")
              ),
              column(
                width = 6,
                numericInput(
                  "pvalue",
                  "Set p.value cutoff",
                  value = rv$params$ora.pv,
                  min = 0, max = 1,
                  step = 0.01
                ),
                shinyBS::bsTooltip("pvalue", 
                          "Maximum threshold for defining subset of genes for ORA analysis", 
                          placement = "bottom", trigger = "hover")
              )
            )
          )
        })
      } # end if ORA
      # dataset is ready if we get to the end of this function
      rv$ds.status <- TRUE
    } # end if gene column
  }
  
  # Stash params
  observeEvent(input$fromtype, {
    rv$params$fromType <- input$fromtype
  })
  observeEvent(input$orgname, {
    rv$params$org.db.name <- input$orgname
  })
  observeEvent(input$minGSSize, {
    rv$params$minGSSize <- input$minGSSize
  })
  observeEvent(input$maxGSSize, {
    rv$params$maxGSSize <- input$maxGSSize
  })
  observeEvent(input$foldchange, {
    rv$params$ora.fc <- input$foldchange
  })
  observeEvent(input$pvalue, {
    rv$params$ora.pv <- input$pvalue
  })
  
  #render table data
  output$table.sample.ds <- DT::renderDataTable(server = TRUE,{
    DT::datatable({ head(getDatasets()) },
                  rownames=FALSE,
                  options = list(
                    dom = 'Brt'
                  )
    )
  })
  
  #update analysis plan
  output$analysis.plan <- renderText(rv$plan)

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
  
  ############
  # analysis #
  ############
  
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
      output$run.header <- renderText("<h3>&nbsp;&nbsp;3. Run Analyses</h3>")
      output$run.button <- renderUI({
        actionButton("run", "Let's go!")
      })
    } else {
      output$run.header <- NULL
      output$run.button <- NULL
      shinyjs::js$hidebox('progress')
    }
  })
  
  # Create progress box
  output$progress <- renderUI({
    box(
        title = "Progress", status = "primary", 
        solidHeader = TRUE, collapsible = T,
        htmlOutput("run_progress"),
        htmlOutput("run.status"),
        uiOutput("view.results")
    )
  })

  ## Button action: run
  observeEvent(input$run, {
    #update UI
    shinyjs::disable("run")
    shinyjs::js$togglebox("db")
    shinyjs::js$togglebox("ds")
    shinyjs::js$showbox('progress')

    # Track Progress 
    withProgress(message = 'Running analyses', value = 0, {
      run.db.list <- strsplit(rv$params$db.list, ",")[[1]]
      run.ds.list <- input$datasets
      steps <- 1 + #libs
        length(run.ds.list) + #proc_datasets
        length(run.ds.list) * length(run.db.list) #run_ora
      i <- 0
      if (rv$params$run.gsea)
        steps <- steps + length(run.ds.list) * length(run.db.list) #run_gsea
      
      # Start time
      start.time <- Sys.time()
      
      # Append progress box html
      prog <- paste("Analysis start time: ",
                    format(start.time, "%F %H:%M:%S"),
                    "<br />",
                    "Steps: ", steps, "<br />")
      rv$logfile <- append(rv$logfile, prog)
      shinyjs::html(id = 'run_progress', add = TRUE, html = prog)
      
      #load libs
      options(install.packages.check.source = "no")
      options(install.packages.compile.from.source = "never")
      if (!require("pacman")) install.packages("pacman"); library(pacman)
      load.libs <- c(load.libs, rv$params$org.db.name)
      p_load(load.libs, update = TRUE, try.bioconductor=TRUE, character.only = TRUE)
      status <- sapply(load.libs,require,character.only = TRUE)
      if(all(status)){
        i <- i + 1
        setProgress(i/steps, detail = paste("Step",i,"of",steps))
        prog <- paste(paste0(i,"."),
          "All required libraries successfully installed and loaded.",
          "<br />")
        rv$logfile <- append(rv$logfile, prog)
        shinyjs::html(id = 'run_progress', add = TRUE, html = prog)
        rv$lib.status <- TRUE
      } else{   
      i <- i + 1
      setProgress(i/steps, detail = paste("Step",i,"of",steps))
        prog <- paste(paste0(i,"."),
          "<p><font color=\"#FF0000\"><b>",
          "ERROR: One or more libraries failed to install correctly.",
          "</b></font></p>",
          "<p>Check your console for FALSE cases and try again...",
          "</p><br /><br />")
        rv$logfile <- append(rv$logfile, prog)
        shinyjs::html(id = 'run_progress', add = TRUE,html = prog)
        rv$lib.status <- FALSE
      }
      # Proceed if libs installed successfully
      if(rv$lib.status){
        #source functions
        source("../scripts/proc_dataset.R")
        source("../scripts/run_ora.R")
        source("../scripts/run_gsea.R")
        source("../scripts/plot_results.R")
        #make output dirs
        output.name <- format(start.time, "%Y%m%d_%H%M%S") # timestamp
        if (dir.exists(file.path("../../shiny_result/")))
          output.name <- file.path("../shiny_result/output",output.name)
        #process each dataset: check errors, generate gene lists, create dirs, etc
        for (ds.name in run.ds.list){
          i <- i + 1
          setProgress(i/steps, detail = paste("Step",i,"of",steps))
          prog <- paste(paste0(i,"."),
            "Preprocessing:",
            ds.name,
            "<br />")
          rv$logfile <- append(rv$logfile, prog)
          shinyjs::html(id = 'run_progress', add = TRUE, html = prog)
          proc_dataset(ds.name, rv$params$org.db.name, rv$params$fromType, 
                       rv$params$ora.fc, rv$params$ora.pv, 
                       rv$params$run.ora, rv$params$run.gsea, output.name)
          #also save params
          ds.noext <- strsplit(ds.name,"\\.")[[1]][1]
          output.dir <- file.path("../",output.name, ds.noext)
          this.fn <- paste0(ds.noext, "__ora_params.rds")
          saveRDS(rv$params, file.path(output.dir,"ora",this.fn)) 
          if(rv$params$run.gsea){
            this.fn <- paste0(ds.noext, "__gsea_params.rds")
            saveRDS(rv$params, file.path(output.dir,"gsea",this.fn)) 
          }

        }
        #run analyses on each dataset
        for (ds.name in run.ds.list){
          for (db.name in run.db.list){
            i <- i + 1
            setProgress(i/steps, detail = paste("Step",i,"of",steps))
            prog <- paste(paste0(i,"."),
              "Running ORA on",
              ds.name, "against", db.name,
              "<br />")
            rv$logfile <- append(rv$logfile, prog)
            shinyjs::html(id = 'run_progress', add = TRUE, html = prog)
            run_ora(ds.name, db.name, output.name)
            if(rv$params$run.gsea){
              i <- i + 1
              setProgress(i/steps, detail = paste("Step",i,"of",steps))
              prog <- paste(paste0(i,"."),
                "Running GSEA on",
                ds.name, "against", db.name,
                "<br />")
              rv$logfile <- append(rv$logfile, prog)
              shinyjs::html(id = 'run_progress', add = TRUE, html = prog)
              run_gsea(ds.name, db.name, output.name)
            }
          }
        }
        # when event complete...
        prog <- paste(
          "<b>Analysis successfully completed:",
          format(Sys.time(), "%F %H:%M:%S"),
          "</b>")
        rv$logfile <- append(rv$logfile, prog)
        shinyjs::html(id = 'run_progress', add = TRUE, html = prog)
        shinyjs::html(id = 'run_progress', add = TRUE, html = paste0(
          "<br /><br /><p>Results can be found in <b>",
          output.name,"/.</b></p>"
        ))
        output$view.results <- renderUI ({
          tagList(
            fluidRow(
              column(width = 4,
                     actionButton("results","View Results")
              ),
              column(width = 4,
                     actionButton("reset","Start new analysis")
              )
            )
          )
        })
        output$run.status <- renderText("&emsp;&#10004; Finished!")
        write(rv$logfile, file.path("../",output.name,"logfile.html"))
      }
    })
  })
  
  observeEvent(input$results, {
    rstudioapi::jobRunScript(path = "../scripts/view_results.R")
  })
  
  observeEvent(input$reset, {
    updateSelectInput(session, "datasets", selected = "")
    shinyjs::js$togglebox("db")
    shinyjs::js$togglebox("ds")
    shinyjs::js$startnew()
  })
  
})
