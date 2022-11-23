#load required packages
require(EnhancedVolcano)
require(writexl)
require(clusterProfiler)
require(ggplot2)
require(enrichplot)

shinyServer(function(input, output, session) {
  
  ########### 
  # sidebar #
  ###########
  
  #update selection lists for drop downs
  observeEvent(input$run,{
    ds.list <- unique(output.df[which(output.df$run==input$run),'dataset'])
    updateSelectizeInput(
      session,
      'dataset', 
      choices = ds.list, 
      server = TRUE
    )
  })
  
  observeEvent(input$dataset,{
    method.list <- unique(output.df[which(output.df$dataset==input$dataset),'method'])
    updateSelectizeInput(
      session,
      'method', 
      choices = method.list, 
      server = TRUE
    )
  })
  
  ############
  # data tab #
  ############ 
  
  #get data params
  getDataParams <- function(){
    fn <- paste(input$dataset, input$method, sep = "__")
    fn <- paste0(fn, "_params.rds")
    fp <- file.path("../output",input$run, input$dataset, input$method, fn)
    readRDS(fp)
  }
  
  #make prefix data filenames (read and download)
  makeDataPrefix <- function(){
    fn <- paste(input$dataset, input$method, sep = "__")
    if(input$excluded.data)
      fn <- paste0(fn, "_excluded")
    else
      fn <- paste0(fn, "_input")
    return(fn)
  }
  
  #get table data
  getTableData <- function(){
    fn <- paste0(makeDataPrefix(), ".rds")
    fp <- file.path("../output",input$run, input$dataset, input$method, fn)
    readRDS(fp)
  }
  
  
  #render table data
  output$table.data <- DT::renderDataTable(server = FALSE,{
    DT::datatable({ getTableData() },
                  rownames=FALSE,
                  extensions = 'Buttons',
                  options = list(
                    dom = 'Bfrtip',
                    pageLength = 20,
                    buttons = list(list(
                      extend = 'collection',
                      buttons = list(
                        list(extend = "excel", title = makeDataPrefix()),
                        list(extend = "csv", title = makeDataPrefix()),
                        list(extend = "copy", title = makeDataPrefix())),
                      text = 'Download'
                    ))
                  )
    )
  })
  
  #plot data
  makePlotData <- function(){
    params <- getDataParams()
    data <- getTableData()
    if('p.value' %in% names(data) & 'fold.change' %in% names(data)) {
      if(!params$fromType %in% names(data)) # i.e., excluded cases
        data[params$fromType] <- data$gene
      selectLab.list <- head(data[[params$fromType]],10)
      if(input$method=='gsea')
        selectLab.list <- c(head(data[[params$fromType]]),
                            tail(data[[params$fromType]]))
      EnhancedVolcano(data,
                      lab = data[[params$fromType]],
                      selectLab = selectLab.list,
                      drawConnectors = TRUE,
                      widthConnectors = 0.2,
                      x = 'fold.change',
                      y = 'p.value',
                      pCutoff = params$pv,
                      FCcutoff = params$fc,
                      legendLabels=c('NS','FC','p-value',
                                     'p-value & FC'),
                      pointSize = 2.0,
                      labSize = 5.0)
    }
  }
  

  #render and cache plot data
  output$plot.data <- renderPlot({
    # ggsave(filename = 'cached.pdf',
    #        plot = makePlotData(),
    #        device = 'pdf',
    #        width = 2400,
    #        height = 2400,
    #        units = "px")
    makePlotData()
    })
  
  #download plot data
  output$download.plot.data <- downloadHandler(
    filename =  function() {
      paste(makeDataPrefix(), "volcano.pdf", sep = "_")
    },
    content = function(file) {
      pdf(file)
      p<-makePlotData()
      print(p)
      dev.off()
    },
    contentType = "application/pdf"
    # filename =  function() {
    #   'cached.pdf'
    # },
    # content = function(file) {
    #   file.copy('cached.pdf',
    #             file, overwrite=TRUE)
    # }
  )

  ##############
  # result tab #
  ##############
  #make prefix data file names (read and download)
  makeResultPrefix <- function(){
    fn <- paste(input$dataset, input$database, input$method, sep = "_")
    return(fn)
  }
  
  # Result Object
  getResultObj <- function(){
    fn <- paste0(makeResultPrefix(), "_result.rds")
    fp <- file.path("../output",input$run, input$dataset, input$method, fn)
    result.obj <- readRDS(fp)
    return(result.obj)
  }
  
  # Gene List
  getGeneList <- function(){
    fn <- paste0(makeResultPrefix(), "_geneList.rds")
    fp <- file.path("../output",input$run, input$dataset, input$method, fn)
    readRDS(fp)
  }
  
  output$debug.text<-renderText({
    # getResultObj()
    # sample.result.id
    # result.obj<-getResultObj()
    # result.obj@result$ID[1]
    # input$table.result_rows_selected
    })
  
  # add database-specific choices
  appendPerDatabase <- function(choices){
    sample.result.id <- getResultObj()@result$ID[1]
    if (grepl("^WP\\d+$", sample.result.id)) #WikiPathways
      unique(c("WikiPathways",choices))
    else if (grepl("^PMC\\d+__", sample.result.id))  #PFOCR
      unique(c("Pathway Figure",choices))
    else if (grepl("^GO:\\d+", sample.result.id))  #GO
      unique(c("Gene Ontology",choices))
  }
  
  #update database-dependent plot options
  observeEvent(input$database, {
    #plot2
    plot2.gsea.choices = c("GSEA score","STRING network")
    plot2.ora.choices = c("STRING network")
    if(input$method == "gsea"){
      plot2.gsea.choices <- appendPerDatabase(plot2.gsea.choices)
      updateSelectInput(session, "plot2", choices = plot2.gsea.choices )
    } else if(input$method == "ora") {
      plot2.ora.choices <- appendPerDatabase(plot2.ora.choices)
      updateSelectInput(session, "plot2", choices = plot2.ora.choices )
    }
  })
  #update method-dependent plot options
  observeEvent(input$method, {
    #plot1
    plot1.ora.choices = c("Dot plot (gene ratio)",
                "Dot plot (count)",
                "Emap plot",
                "Concept network",
                "Concept network (circular)",
                "Volcano plot (ORA)",
                "Heatmap (ORA)",
                "Upset (ORA)")
    plot1.gsea.choices = c("Dot plot (gene ratio)",
                          "Dot plot (count)",
                          "Emap plot",
                          "Concept network",
                          "Concept network (circular)",
                          "Volcano plot (GSEA)",
                          "Heatmap (GSEA)")
    plot1.other.choices = c("Dot plot (gene ratio)",
                           "Dot plot (count)",
                           "Emap plot",
                           "Concept network",
                           "Concept network (circular)")
    if(input$method == "gsea")
      updateSelectInput(session, "plot1", choices = plot1.gsea.choices )
    else if(input$method == "ora")
      updateSelectInput(session, "plot1", choices = plot1.ora.choices )
    else
      updateSelectInput(session, "plot1", choices = plot1.other.choices )
 
   #plot2
    plot2.gsea.choices = c("GSEA score","STRING network")
    plot2.ora.choices = c("STRING network")
    if(input$method == "gsea"){
      plot2.gsea.choices <- appendPerDatabase(plot2.gsea.choices)
      updateSelectInput(session, "plot2", choices = plot2.gsea.choices )
    } else if(input$method == "ora") {
      plot2.ora.choices <- appendPerDatabase(plot2.ora.choices)
      updateSelectInput(session, "plot2", choices = plot2.ora.choices )
    }
  })
  
  #render table result
  output$table.result <- DT::renderDataTable(
    server = FALSE, # FALSE to support download of all table rows
    withProgress(message = 'Updating table...', detail = "",
                 style = getShinyOption("progress.style", default = "notification"),
                 value = NULL, { 
    DT::datatable({ as.data.frame(getResultObj()) },
                  rownames=FALSE,
                  selection = list(mode = 'single', selected = c(1)),
                  extensions = 'Buttons',
                  options = list(
                    dom = 'Bfrtip',
                    pageLength = 10,
                    # scrollX=TRUE,
                    search = list(regex = TRUE, caseInsensitive = TRUE),
                    columnDefs = list(
                      list(targets = c(7:10), visible = FALSE, defaultContent = "-"),
                      list(targets=c(1), visible=TRUE, width='250',  #description
                           render = JS(
                             "function(data, type, row, meta) {",
                             "return type === 'display' && data != null && data.length > 50 ?",
                             "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
                             "}")
                      )
                    ),
                    buttons = list(list(
                      extend = 'collection',
                      buttons = list(
                        list(extend = "excel", title = makeResultPrefix()),
                        list(extend = "csv", title = makeResultPrefix()),
                        list(extend = "copy", title = makeResultPrefix())),
                      text = 'Download'
                    ))
                  )
    )
   })
  )
  
  #plot1 result
  makePlot1Result <- function(){
    resObject <- getResultObj()
    data.volcano <- as.data.frame(resObject)
    data.emap <- pairwise_termsim(resObject)
    switch (input$plot1,
            "Dot plot (gene ratio)" = dotplot(resObject,
                                         showCategory = 20,
                                         label_format=50),
            "Dot plot (count)" = dotplot(resObject,
                                    showCategory = 20,
                                    x = "count",
                                    label_format=50),
            "Emap plot" = emapplot(data.emap, 
                                   showCategory = 20,
                                   cex_label_category=0.7,
                                   layout="nicely"), #alt layouts: "kk","sugiyama","nicely","fr", "gem","lgl","mds",
            "Concept network" = cnetplot(resObject, 
                                         foldChange=getGeneList(),
                                         categorySize="geneNum", 
                                         cex_label_category = 0.8, 
                                         cex_label_gene = 1.0),
            "Concept network (circular)" = cnetplot(resObject, 
                                                    foldChange=getGeneList(), 
                                                    circular = TRUE, 
                                                    colorEdge = TRUE, 
                                                    cex_label_category = 0.8, 
                                                    cex_label_gene = 1.0) ,
            "Volcano plot (GSEA)" = EnhancedVolcano(data.volcano,
                                             lab = data.volcano$Description,
                                             selectLab = head(data.volcano$Description,3),
                                             drawConnectors = TRUE,
                                             widthConnectors = 0.2,
                                             x = 'NES',
                                             y = 'pvalue',
                                             pCutoff = 0.05,
                                             FCcutoff = 1,
                                             legendLabels=c('NS','NES','p-value',
                                                            'p-value & NES'),
                                             xlab = 'Normalized Enrichment Score',
                                             pointSize = 2.0,
                                             labSize = 5.0),
            "Volcano plot (ORA)" = EnhancedVolcano(data.volcano,
                                                   lab = data.volcano$Description,
                                                   selectLab = head(data.volcano$Description,3),
                                                   drawConnectors = TRUE,
                                                   widthConnectors = 0.2,
                                                   x = 'Count',
                                                   y = 'pvalue',
                                                   pCutoff = 0.05,
                                                   FCcutoff = 5,
                                                   legendLabels=c('NS','Count','p-value',
                                                                  'p-value & Count'),
                                                   xlab = 'Number of Genes',
                                                   pointSize = 2.0,
                                                   labSize = 5.0),
            "Heatmap (GSEA)" = heatplot(resObject, foldChange=getGeneList(), 
                                 showCategory=5,
                                 label_format=50) + coord_fixed(ratio=2),
            "Heatmap (ORA)" = heatplot(resObject, #can't figure out foldChange param
                                       showCategory=5,
                                       label_format=50) + coord_fixed(ratio=2),
            "Upset plot (ORA)" = upsetplot(resObject, n=10)
    )
  }
  
  
  
  #render and cache plot data
  output$plot1.result <- renderPlot({
    # ggsave(filename = 'cached1.pdf',
    #        plot = makePlot1Result(),
    #        device = 'pdf',
    #        width = 2400,
    #        height = 2400,
    #        units = "px")
    makePlot1Result()
  })
  
  #download plot data
  output$download.plot1.result <- downloadHandler(
    filename =  function() {
      paste0(paste(makeResultPrefix(),
            gsub("\\(|\\)","",gsub(" ","-",input$plot1)), 
            sep = "_"), ".pdf")
    },
    content = function(file) {
      pdf(file)
      p<-makePlot1Result()
      print(p)
      dev.off()
    },
    contentType = "application/pdf"
    # filename =  function() {
    #   'cached1.pdf'
    # },
    # content = function(file) {
    #   file.copy('cached1.pdf',
    #             file, overwrite=TRUE)
    # }
  )
  
  #html result
  makeHtmlResult <- function(){
    resObject <- getResultObj()
    req(input$table.result_rows_selected) #wait for table to load
    i <- input$table.result_rows_selected #row selection
    res.object.id <- resObject$ID[i]
    res.object.id <- gsub("\\.jpg$","",res.object.id) #pfocr ids
    switch (input$plot2,
            "GSEA score" = NULL,
            "STRING network" = NULL,
            "WikiPathways" = paste0('View <a href="https://new.wikipathways.org/pathways/',
                                           res.object.id,
                                           '" target="_blank">',
                                           res.object.id,
                                           ' at WikiPathways</a><br/><br/>'),
            "Pathway Figure" = paste0('View <a href="https://pfocr.wikipathways.org/figures/',
                                             res.object.id,
                                             '" target="_blank">',
                                             res.object.id,
                                             ' at PFOCR</a><br/><br/>'),
            "Gene Ontology" = paste0('View <a href="http://amigo.geneontology.org/amigo/term/',
                                            res.object.id,
                                            '" target="_blank">',
                                            res.object.id,
                                            ' at Gene Ontology</a><br/><br/>')
    )
  }
  
  #render html
  output$html.result <- renderText({ 
    makeHtmlResult() 
    })
  
  #plot2 result
  makePlot2Result <- function(){
    resObject <- getResultObj()
    req(input$table.result_rows_selected) #wait for table to load
    i <- input$table.result_rows_selected #row selection
    res.object.id <<- resObject$ID[i]
    switch (input$plot2,
            "GSEA score" = gseaplot2(resObject, 
                                     geneSetID = i, 
                                     title = resObject$Description[i]) ,
            "STRING network" = NULL,
            "WikiPathways" = NULL,
            "Pathway Figure" = NULL,
            "Gene Ontology" = NULL
    )
  }
  
  #render and plot 
  output$plot2.result <- renderPlot({
    makePlot2Result()
  })
  
  #download plot data
  output$download.plot2.result <- downloadHandler(
    filename =  function() {
      paste0(paste(makeResultPrefix(),
            gsub("\\(|\\)","",gsub(" ","-",input$plot2)), 
            gsub(" ","-",res.object.id), 
            sep = "_"), ".pdf")
    },
    content = function(file) {
      pdf(file)
      p<-makePlot2Result()
      print(p)
      dev.off()
    },
    contentType = "application/pdf"
  )

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
