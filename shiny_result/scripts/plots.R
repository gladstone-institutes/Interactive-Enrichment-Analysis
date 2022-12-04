# All the plotting functions used in shiny_result
## Add or edit individual plots and their options in PLOTS & OPTIONS
## Add or edit the collection of individual option inputs in OPTION INPUTS

library(shinyjs)
library(stringr)
library(dplyr)
library(ggplot2)
library(enrichplot)

################################################################################
# OPTION INPUTS #
#################
# Reusable option inputs are all named here (all.options) and the filterOptions
# function is called from each PLOTS & OPTIONS function to specify the subset
# to be shown.

filterOptions <- function(options.list){
  all.options <- c("showCategory","showGene","label_format","dotplot_x",
                   "dotplot_size","category_label","category_node",
                   "gene_label","gene_node","color","layout", "cnet_layout",
                   "cnet_circular", "cnet_colorEdge","bar_by")
  sapply(all.options, shinyjs::show)
  sapply(setdiff(all.options,options.list), shinyjs::hide)
}

################################################################################
# PLOTS & OPTIONS #
###################
# Each function returns a plot object for rendering or printing/saving and 
# defines an option panel with corresponding inputs.

## Pro-tip: 


# Dotplot #
# - standard enrichplot function with "GeneRatio" or "Count" x-axis option
shinyDotplot <- function(resObject, input, output){
  
  #filter inputs options panel
  filterOptions(c("showCategory","label_format","dotplot_x","dotplot_size",
                  "category_label","color"))
  #initalize
  size = input$dotplot_size
  
  #Apparently not implemented by enrichplot...
  if(size == "Percentage"){
    if(input$method == "ora"){
      resObject@result <- resObject@result %>%
        dplyr::rowwise() %>%
        dplyr::mutate(Percentage = Count/as.integer(
          str_split(BgRatio, "\\/")[[1]][1]))
    } else if (input$method == "gsea"){
      resObject@result <- resObject@result %>%
        dplyr::rowwise() %>%
        dplyr::mutate(Percentage = stringr::str_count(core_enrichment, "\\/")/setSize)
    }
  }
  
  #plot
  enrichplot::dotplot(resObject,
                      showCategory = input$showCategory,
                      label_format=input$label_format,
                      font.size = input$category_label,
                      size = size,
                      color = input$color,
                      x = input$dotplot_x)
}

# Emapplot #
# - standard enrichplot function with layout options
shinyEmapplot <- function(data.emap, input, output){
  
  #filter inputs options panel
  filterOptions(c("showCategory","label_format","category_label",
                  "category_node","color","layout"))
  
  #plot
  enrichplot::emapplot(data.emap, 
                       showCategory = input$showCategory,
                       color = input$color,
                       layout.params=list(layout=input$layout),
                       cex.params = list(
                         category_label=input$category_label/12, #pt to rem
                         category_node = input$category_node
                       )
  )
}

# Cnetplot #
# - standard enrichplot function with circular option
shinyCnetplot <- function(resObject, geneList, input, output){
  
  #filter inputs options panel
  filterOptions(c("showCategory","category_label","gene_label","category_node",
                  "gene_node","cnet_layout", "cnet_circular", "cnet_colorEdge"))
  
  #plot
  enrichplot::cnetplot(resObject, 
                       showCategory = input$showCategory,
                       circular = input$cnet_circular, 
                       colorEdge = input$cnet_colorEdge, 
                       layout = input$cnet_layout,
                       cex.params = list(
                         foldChange=geneList,
                         category_label=input$category_label/12, #pt to rem
                         category_node = input$category_node,
                         gene_label = input$gene_label/12,
                         gene_node = input$gene_node
                       )
  )
}

# Heatmap #
# - starts with top n results, identifies the top n genes by frequency within
# those results (secondary sort on fold.change) and then makes a tile plot with 
# fold.change fill color using a balanced Brewer RdBu pallette
shinyHeatmap <- function(resObject, data, params, input, output){
  
  #filter inputs options panel
  filterOptions(c("showCategory","showGene","category_label","gene_label"))
  
  #prep data
  res <- resObject@result[1:input$showCategory, c('Description','geneID')]
  res <- tidyr::separate_rows(res, geneID, sep="\\/")
  #define secondary order and color mapping if available
  if('fold.change' %in% names(data)) {
    res <- dplyr::left_join(res,data, by=c("geneID"=params$fromType))
  } else {
    res$fold.change = 0
  }
  res <- dplyr::add_count(res, geneID)
  res <- dplyr::arrange(res, desc(n),desc(abs(fold.change)))
  freq.gene <- unique(res$geneID)
  res <- dplyr::filter(res, geneID %in% freq.gene[1:input$showGene]) 
  res <- dplyr::mutate(res, Description = stringr::str_trunc(Description, 50))
  res <- dplyr::mutate(res, geneID = factor(
    geneID, levels = unique(geneID))) #fix two-factor sorting
  res <- dplyr::mutate(res, Description = factor(
    Description, levels=unique(Description)))
  
  #plot
  p <- ggplot(res, aes(x=geneID, 
                       y=Description, fill = fold.change)) + 
    geom_tile(color = 'white') +
    xlab(NULL) + ylab(NULL) + theme_minimal() +
    theme(panel.grid.major = element_blank(),
          legend.position = "none",
          axis.text.x=element_text(angle = 60, hjust = 1, 
                                   size=input$gene_label),
          axis.text.y=element_text(size=input$category_label))
  if('fold.change' %in% names(data)){
    max.scale <- max(abs(res$fold.change))
    p <- p + scale_fill_distiller(type = "div",
                                  palette = "RdBu",
                                  direction = -1, 
                                  limits = c(-max.scale, max.scale)) +
      theme(legend.position = "right")
  }
  return(p)
}

# BarplotResult #
# - ...
shinyBarplotResult <- function(resObject, data, params, input, output){
  
  #filter inputs options panel
  filterOptions(c("showCategory","category_label","bar_by"))
  
  #prep data 
  res <- resObject@result[1:input$showCategory, ]
  res <- tidyr::separate_rows(res, geneID, sep="\\/")
  #define secondary order and color mapping if available
  if('fold.change' %in% names(data)) {
    res <- dplyr::left_join(res,data, by=c("geneID"=params$fromType))
  } else {
    res$fold.change = 0
  }
  if(input$method == "ora"){
    res <- res %>%
      dplyr::group_by(Description, BgRatio, Count) %>%
      dplyr::summarise("Fold change" = mean(fold.change)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Percentage = Count/as.integer(
        str_split(BgRatio, "\\/")[[1]][1]))
  } else if (input$method == "gsea"){
    res <- res %>%
      dplyr::add_count(Description, name = "Count") %>%
      dplyr::group_by(Description, Count, setSize) %>%
      dplyr::summarise("Fold change" = mean(fold.change)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Percentage = Count/setSize)
  }
  res<- res %>%
    dplyr::arrange( !!as.name(input$bar_by)) %>%
    dplyr::mutate(Description = factor(
      Description, levels=unique(Description)))
  
  #create a column with positive/negative expressed genes
  res$DEG <- NA
  if(input$bar_by == 'Fold change'){
    res$DEG[res$`Fold change` > 0] <- "upregulated"
    res$DEG[res$`Fold change` < 0] <- "downregulated"
  }
  
  #plot
  ggplot(res, 
         aes(!!as.name(input$bar_by), Description, fill=DEG)) +
    geom_bar(stat="identity") +
    # ggbreak::scale_y_break(c( -3, -5.9), scale=3)+ 
    scale_fill_manual(values=c("#67A9CF","#EF8A62"),
                      name = element_blank()) +
    guides(fill = guide_legend(reverse = TRUE)) +
    xlab(input$bar_by) + ylab("") +
    theme_bw() +
    theme(text = element_text(size = input$category_label),
          axis.title.y = element_blank(),
          legend.position = "right")
  
  # res <- dplyr::add_count(res, geneID)
  # res <- dplyr::arrange(res, desc(n),desc(abs(fold.change)))
  # res <- dplyr::mutate(res, Description = stringr::str_trunc(Description, 50))
  # res <- dplyr::mutate(res, geneID = factor(
  #   geneID, levels = unique(geneID))) #fix two-factor sorting
  # res <- dplyr::mutate(res, Description = factor(
  #   Description, levels=unique(Description)))
  # 
 
  
  #plot
  # p <- ggplot(res, aes(x=geneID, 
  #                      y=Description, fill = fold.change)) + 
  #   geom_tile(color = 'white') +
  #   xlab(NULL) + ylab(NULL) + theme_minimal() +
  #   theme(panel.grid.major = element_blank(),
  #         legend.position = "none",
  #         axis.text.x=element_text(angle = 60, hjust = 1, 
  #                                  size=input$gene_label),
  #         axis.text.y=element_text(size=input$category_label))
  # if('fold.change' %in% names(data)){
  #   max.scale <- max(abs(res$fold.change))
  #   p <- p + scale_fill_distiller(type = "div",
  #                                 palette = "RdBu",
  #                                 direction = -1, 
  #                                 limits = c(-max.scale, max.scale)) +
  #     theme(legend.position = "right")
  # }
  # return(p)
}

##############
# DATA PLOTS #
##############

# Volcano #
# - an EnhancedVolcano plot of pv and fc, exposing options
# for top n genes (or gene list), font size and legend position.
shinyVolcano <- function(data, params, input, output){
  
  #handle selection
  selectLab.list <- input$selectLab
  if(length(selectLab.list) == 0){
    selectLab.list <- head(data[[params$fromType]], input$topLab)
    if(input$method=='gsea')
      selectLab.list <- c(head(data[[params$fromType]],as.integer(input$topLab/2)),
                          tail(data[[params$fromType]],as.integer(input$topLab/2)))
  } 
  
  #plot 
  EnhancedVolcano(data,
                  lab = data[[params$fromType]],
                  selectLab = selectLab.list,
                  drawConnectors = TRUE,
                  widthConnectors = 0.2,
                  x = 'fold.change',
                  y = 'p.value',
                  title = "",
                  subtitle = "",
                  legendPosition = input$legend_pos,
                  pCutoff = params$ora.pv,
                  FCcutoff = params$ora.fc,
                  legendLabels=c('NS','FC','p-value',
                                 'p-value & FC'),
                  pointSize = 2.0,
                  labSize = input$category_label_0/3)
}


# Barplot #
# - standard bar plot distinguishing up/down regulated genes, exposing options
# for top n genes (or gene list), font size and legend position.
shinyBarplot <- function(data, params, input, output){
  
  #handle selection
  selectLab.list <- input$selectLab
  if(length(selectLab.list) == 0){
    selectLab.list <- head(data[[params$fromType]], input$topLab)
    if(input$method=='gsea')
      selectLab.list <- c(head(data[[params$fromType]],as.integer(input$topLab/2)),
                          tail(data[[params$fromType]],as.integer(input$topLab/2)))
  } 
  #subset data accordingly
  data <- dplyr::filter(data, !!as.name(params$fromType) %in% selectLab.list)
  
  #create a column with positive/negative expressed genes
  data$DEG <- NA
  data$DEG[data$fold.change>0] <- "upregulated"
  data$DEG[data$fold.change<0] <- "downregulated"
  
  #plot 
  ggplot(data, 
         aes(selectLab.list, fold.change, fill=DEG)) +
    geom_bar(stat="identity") +
    # ggbreak::scale_y_break(c( -3, -5.9), scale=3)+ 
    scale_fill_manual(values=c("#67A9CF","#EF8A62"),
                      name = element_blank()) +
    guides(fill = guide_legend(reverse = TRUE)) +
    xlab("Genes") + ylab(expression(Log[2]~fold~change)) +
    theme_bw() +
    theme(text = element_text(size = input$category_label_0),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = input$legend_pos)
}

