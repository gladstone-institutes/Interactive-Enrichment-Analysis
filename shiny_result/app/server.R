library(shinyjs)
library(stringr)
library(writexl)
library(ggplot2)
library(enrichplot)
library(EnhancedVolcano)
source("../scripts/plots.R", local = TRUE)


shinyServer(function(input, output, session) {
  
  ########### 
  # sidebar #
  ###########
  
  #update selection lists for drop downs
  observeEvent(input$run,{
    ds.list <- unique(output.df[which(output.df$run==input$run),'dataset'])
    method.list <- rev(unique(output.df[which(output.df$run==input$run &
                                            output.df$dataset==input$dataset),'method']))
    updateSelectizeInput(
      session,
      'dataset',
      choices = ds.list,
      server = TRUE
    )
    updateSelectizeInput(
      session,
      'method',
      choices = method.list,
      server = TRUE
    )
  })
  
  observeEvent(input$dataset,{
    method.list <- rev(unique(output.df[which(output.df$run==input$run &
                                            output.df$dataset==input$dataset),'method']))
    # db.list <- strsplit(getDataParams()$db.list, ",")[[1]]
    updateSelectizeInput(
      session,
      'method',
      choices = method.list,
      server = TRUE
    )
    updateSelectizeInput(
      session,
      'databse',
      choices = db.list,
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
                  selection = "none",
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
        selectLab.list <- c(head(data[[params$fromType]],5),
                            tail(data[[params$fromType]]),5)
      #create a column with positive/negative expressed genes
      data$DEG <- NA
      data$DEG[data$fold.change>0] <- "upregulated"
      data$DEG[data$fold.change<0] <- "downregulated"
      switch (input$plot0,
              "Volcano plot" = EnhancedVolcano(data,
                      lab = data[[params$fromType]],
                      selectLab = selectLab.list,
                      drawConnectors = TRUE,
                      widthConnectors = 0.2,
                      x = 'fold.change',
                      y = 'p.value',
                      pCutoff = params$ora.pv,
                      FCcutoff = params$ora.fc,
                      legendLabels=c('NS','FC','p-value',
                                     'p-value & FC'),
                      pointSize = 2.0,
                      labSize = 5.0),
              "Heatmap" = NULL,
              "Bar plot" = ggplot(head(data,10), 
                                  aes(head(data[[params$fromType]],10), fold.change, fill=DEG)) +
                geom_bar(stat="identity") +
                # ggbreak::scale_y_break(c( -3, -5.9), scale=3)+ 
                scale_fill_manual(values=c("#67A9CF","#EF8A62"),
                                  name = element_blank()) +
                guides(fill = guide_legend(reverse = TRUE)) +
                xlab("Genes") + ylab(expression(Log[2]~fold~change)) +
                theme_bw() +
                theme(text = element_text(size = 16),
                      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
      )
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
    gl <- readRDS(fp)
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
    plot2.gsea.choices = c("GSEA score","Linkouts","STRING network")
    plot2.ora.choices = c("Linkouts","STRING network")
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
    plot1.ora.choices = c("Dot plot",
                "Emap plot",
                "Concept network",
                "Heatmap",
                "Upset (ORA)")
    plot1.gsea.choices = c("Dot plot",
                          "Emap plot",
                          "Concept network",
                          "Heatmap")
    plot1.other.choices = c("Dot plot",
                           "Emap plot",
                           "Concept network")
    if(input$method == "gsea")
      updateSelectInput(session, "plot1", choices = plot1.gsea.choices )
    else if(input$method == "ora")
      updateSelectInput(session, "plot1", choices = plot1.ora.choices )
    else
      updateSelectInput(session, "plot1", choices = plot1.other.choices )
 
   #plot2
    plot2.gsea.choices = c("GSEA score","Linkouts","STRING network")
    plot2.ora.choices = c("Linkouts","STRING network")
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
    #trim descriptions to 80 characters
    resObject@result <- dplyr::mutate(resObject@result, Description = stringr::str_trunc(Description, 80))
    switch (input$plot1,
            "Dot plot" = shinyDotplot(resObject,input, output),
            "Emap plot" = {
              data.emap <- pairwise_termsim(resObject)
              shinyEmapplot(data.emap,input,output)
              },
            "Concept network" = enrichplot::cnetplot(resObject, 
                                         showCategory = input$showCategory,
                                         foldChange=getGeneList(),
                                         categorySize="geneNum", 
                                         cex_label_category = 0.8, 
                                         cex_label_gene = 1.0),
            "Concept network (circular)" = enrichplot::cnetplot(resObject, 
                                                    showCategory = input$showCategory,
                                                    foldChange=getGeneList(), 
                                                    circular = TRUE, 
                                                    colorEdge = TRUE, 
                                                    cex_label_category = 0.8, 
                                                    cex_label_gene = 1.0) ,
            "Heatmap" = {
              data <- getTableData()
              params <- getDataParams()
              if (input$method == "gsea")
                resObject@result$geneID <- resObject@result$core_enrichment 
              shinyHeatmap(resObject, data, params, input, output)
              },
            "Upset plot (ORA)" = enrichplot::upsetplot(resObject, 
                                                       n=input$showCategory)
    )
  }
  
  #render and cache plot data
  output$plot1.result <- renderPlot({
    # ggplot2::ggsave(filename = 'cached1.pdf',
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
      pdf(file, width = input$plot1.width/72, height = input$plot1.height/72) #pixels to inches
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
    #get result ID for selected row
    res.object.id <- resObject$ID[i]
    res.object.id <- gsub("\\.jpg$","",res.object.id) #pfocr ids
    #get gene list for selected row
    res.object.geneList <- ""
    if(input$method == "ora"){
      res.object.geneList <- resObject$geneID[i]
    } else if (input$method == "gsea"){
      res.object.geneList <- resObject$core_enrichment[i]
    }
    res.object.geneList <- strsplit(res.object.geneList, "\\/")[[1]]
    #arbitrary cutoff for GET-based gene list queries
    res.object.geneList.GET <- res.object.geneList
    if(length(res.object.geneList.GET) > 100) 
      res.object.geneList.GET <- res.object.geneList.GET[1:100]
    res.object.geneList.GET.csv <- paste(res.object.geneList, collapse = ",") 
    #define gene lists and color mapping
    data.mapping <- NULL
    params <- getDataParams()
    data <- getTableData()
    #get up/down regulated if available
    if('p.value' %in% names(data) & 'fold.change' %in% names(data)) {
      entrez.up <- data %>%
        dplyr::filter(!!as.name(params$fromType) %in% res.object.geneList) %>%
        dplyr::filter(fold.change > params$ora.fc & p.value < params$ora.pv) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(new.entrez = paste("Entrez","Gene",ENTREZID, sep = "_")) %>% 
        dplyr::pull(new.entrez) #always available in data
      entrez.dn <- data %>%
        dplyr::filter(!!as.name(params$fromType) %in% res.object.geneList) %>%
        dplyr::filter(fold.change < -params$ora.fc & p.value < params$ora.pv) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(new.entrez = paste("Entrez","Gene",ENTREZID, sep = "_")) %>%
        dplyr::pull(new.entrez) #always available in data
      entrez.up <- paste(entrez.up, collapse = ",")
      entrez.dn <- paste(entrez.dn, collapse = ",")
      data.mapping <- paste0('&EF8A62=', entrez.up,
                             '&67A9CF=', entrez.dn)
    } else {
      entrez.all <- data %>%
        dplyr::filter(!!as.name(params$fromType) %in% res.object.geneList.GET) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(new.entrez = paste("Entrez","Gene",ENTREZID, sep = "_")) %>% 
        dplyr::pull(new.entrez) #always available in data
      entrez.all <- paste(entrez.all, collapse = ",")
      data.mapping <- paste0('&lightgreen=', entrez.all)
    }
    #construct custom linkout buttons
    custom.linkout.button <- ""
    if (grepl("^WP\\d+$", res.object.id)) { #WikiPathways
      custom.linkout.button <- makeLinkoutButton(
        "btn-primary",
        "https://upload.wikimedia.org/wikipedia/commons/3/34/Wplogo_500.png",
        "View the selected pathway at WikiPathways",
        paste0('https://new.wikipathways.org/pathways/',
               res.object.id),
        "WikiPathways")
    } else if (grepl("^PMC\\d+__", res.object.id)) { #PFOCR
      custom.linkout.button <- makeLinkoutButton(
        "btn-primary",
        "https://upload.wikimedia.org/wikipedia/commons/3/34/Wplogo_500.png",
        "View the selected pathway figure at PFOCR",
        paste0('https://pfocr.wikipathways.org/figures/',
               res.object.id),
        "Pathway Figures")
    } else if (grepl("^GO:\\d+", res.object.id)) { #GO
      custom.linkout.button <- makeLinkoutButton(
        "btn-primary",
        "https://avatars1.githubusercontent.com/u/7750835?s=200&v=4",
        "View the selected GO term at AmiGO",
        paste0('http://amigo.geneontology.org/amigo/term/',
               res.object.id),
        "Gene Ontology")
    }
    switch (input$plot2,
            "GSEA score" = NULL,
            "STRING network" = NULL,
            "WikiPathways" = paste0(
              custom.linkout.button,
              '<iframe src ="https://pathway-viewer.toolforge.org/?id=',
              res.object.id,
              data.mapping,
              '" height="350px" width="450px" style="overflow:hidden;"></iframe> ',
              '<br /><a href="https://pathway-viewer.toolforge.org/?id=',
              res.object.id,
              data.mapping,
              '" target="_blank" >open in new window</a>'
            ),
            "Linkouts" = paste0(
              '<h3>Linkouts</h3>',
              custom.linkout.button,
              # '<li><a class="drugstone-button drugstone-green" ',
              makeLinkoutButton(
                "btn-primary",
                "https://cdn.drugst.one/libs/drugstone-buttons/0.0.1/android-chrome-192x192.png",
                "Query Drugst.One with the genes from selected pathway",
                paste0('https://drugst.one/standalone?nodes=',
                       res.object.geneList.GET.csv,
                       '&autofillEdges=true&activateNetworkMenuButtonAdjacentDrugs=true&interactionDrugProtein=NeDRex&licensedDatasets=true'),
                "Drugst.One")),
            "Pathway Figure" = paste0(
              custom.linkout.button,
              '<img src="https://www.ncbi.nlm.nih.gov/pmc/articles/',
              strsplit(res.object.id, "__")[[1]][1],
              '/bin/',
              strsplit(res.object.id, "__")[[1]][2],
              '.jpg" style="max-width:450px;">'
            ),
            "Gene Ontology" = paste0(
              custom.linkout.button,
              '<img src="http://amigo.geneontology.org/visualize?term_data_type=string&mode=amigo&inline=false&term_data=',
              res.object.id,
              '&format=png" style="max-width:450px;">',
              '<br /><a href="http://amigo.geneontology.org/visualize?term_data_type=string&mode=amigo&inline=false&term_data=',
              res.object.id,
              '&format=svg" target="_blank" >open in new window</a>'
            )
    )
  }
  
  #custom linkout button
  makeLinkoutButton <- function(btn.class,
                                btn.img.url,
                                btn.tooltip,
                                btn.href,
                                btn.label){
    paste0(
      '<a class="btn ',btn.class,'" ',
      'title="',btn.tooltip,'" ',
      'style="background-image: url(',btn.img.url,');',
      'background-repeat: no-repeat;',
      'background-position: left;',
      'background-position-x: left;',
      'background-size: 20px 20px;',
      'background-position-x: 6px;',
      'background-origin: padding-box;',
      'padding-left: 32px !important; ',
      'margin: 10px;" ',
      'href="',btn.href,
      '" target="_blank">',btn.label,'</a><br />')
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
            "GSEA score" = enrichplot::gseaplot2(resObject, 
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
      grDevices::pdf(file)
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
