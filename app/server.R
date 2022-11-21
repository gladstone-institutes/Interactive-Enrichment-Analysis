#load required packages
require(gridExtra)
require(EnhancedVolcano)
require(writexl)

shinyServer(function(input, output, session) {
  
  ########### 
  # sidebar #
  ###########
  
  #update selection lists for drop downs
  inputgenes <- reactiveValues()
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
                  options = list(dom = 'Bfrtip',
                  pageLength = 20,
                  buttons = list(list(extend = 'collection',
                                      buttons = list(
                                        list(extend = "excel", title = makeDataPrefix()),
                                        list(extend = "csv", title = makeDataPrefix()),
                                        list(extend = "copy", title = makeDataPrefix())),
                                      text = 'Download'
                  ))
                  )
    )
  })
  # output$debug.text<-renderText(getDataParams()$fromType)
  
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
    ggsave(filename = 'cached.pdf',
           plot = makePlotData(),
           device = 'pdf',
           width = 2400,
           height = 2400,
           units = "px")
    makePlotData()
    })
  
  #download plot data
  output$download.plot.data <- downloadHandler(
    filename =  function() {
      'cached.pdf'
    },
    content = function(file) {
      file.copy('cached.pdf',
                file, overwrite=TRUE)
    }
  )
  
  # ggsave(p, file = file.path(output.dir,method.dir,"plots",vol.fn), 
  #        width = 2400, height = 2400, units = "px", device='pdf')

  ##############
  # result tab #
  ##############
  # 
  # #feature plot tab - download button for reference UMAPs
  # refUMAPPlotfunction <- function(){
  #   p1 <- dittoDimPlot(get(input$dataset),
  #                      c("seurat_clusters"),
  #                      order = "increasing",
  #                      size = 2.5,
  #                      opacity = 0.8,
  #                      show.others = T,
  #                      main = "Seurat clusters",
  #                      theme = theme_bw(base_size = 15)
  #                     )
  #   p2 <- dittoDimPlot(get(input$dataset),
  #                      c("cell_type"),
  #                      order = "increasing",
  #                      size = 2.5,
  #                      opacity = 0.8,
  #                      show.others = T,
  #                      main = "Cell type",
  #                      theme = theme_bw(base_size = 15)
  #                     )
  # 
  #   grid.arrange(p1,p2, ncol = 2)
  # }
  # output$downloadFefUMAPPlot <- downloadHandler(
  #   filename =  function() {
  #     #downloaded file name
  #     paste0(plot_download_prefix, "_", Sys.time(),".png")
  #   },
  #   #content is a function with argument file.
  #   #content writes the plot to the device
  #   content = function(file) {
  #     png(file,
  #         width = 580*3,
  #         height = 500)
  #     refUMAPPlotfunction() # draw the plot
  #     dev.off()  # turn the device off
  #   }
  # )
  # 
  # #feature plot tab - reference UMAP plot output (group by seurat cluster)
  # output$refUMAPcluster <- renderPlot({
  #   dittoDimPlot(get(input$dataset),
  #                c("seurat_clusters"),
  #                order = "increasing",
  #                size = 2.5,
  #                opacity = 0.8,
  #                show.others = T,
  #                main = "Seurat clusters",
  #                theme = theme_bw(base_size = 15)
  #   )
  # })
  # 
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
