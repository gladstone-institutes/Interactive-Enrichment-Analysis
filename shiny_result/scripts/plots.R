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
# Reusable option inputs are defined here as functions to be called from 
# one or more plot&options functions, where option panels are composed.

filterOptions <- function(options.list){
  all.options <- c("showCategory","showGene","label_format","dotplot_x",
                   "dotplot_size","cex_label_category","cex_category",
                   "cex_label_gene","color","layout")
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
  filterOptions(c("showCategory","label_format","dotplot_x","dotplot_size","cex_label_category","color"))
  
  #plot
  size = input$dotplot_size
  if(size == "Percentage"){
    resObject@result <- resObject@result %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Percentage = Count/as.integer(str_split(BgRatio, "\\/")[[1]][1]))
  }
  enrichplot::dotplot(resObject,
                      showCategory = input$showCategory,
                      label_format=input$label_format,
                      font.size = input$cex_label_category,
                      size = size,
                      color = input$color,
                      x = input$dotplot_x)
}

# Emapplot #
# - standard enrichplot function with layout options
shinyEmapplot <- function(data.emap, input, output){
  
  #filter inputs options panel
  filterOptions(c("showCategory","label_format","cex_label_category","cex_category","color","layout"))
  
  #plot
  enrichplot::emapplot(data.emap, 
                     showCategory = input$showCategory,
                     cex_label_category=input$cex_label_category/12, #pt to rem
                     cex_category = input$cex_category,
                     color = input$color,
                     layout.params=list(layout=input$layout) 
  )
}

# Heatmap #
# - starts with top n results, identifies the top n genes by frequency within
# those results (secondary sort on fold.change) and then makes a tile plot with 
# fold.change fill color using a balanced Brewer RdBu pallette
shinyHeatmap <- function(resObject, data, params, input, output){
  
  #filter inputs options panel
  filterOptions(c("showCategory","showGene","cex_label_category","cex_label_gene"))
  
  #plot  
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
  res <- dplyr::mutate(res, geneID = factor(geneID, levels = unique(geneID))) #fix two-factor sorting
  res <- dplyr::mutate(res, Description = factor(Description, levels=unique(Description)))
  p <- ggplot(res, aes(x=geneID, 
                       y=Description, fill = fold.change)) + 
    geom_tile(color = 'white') +
    xlab(NULL) + ylab(NULL) + theme_minimal() +
    theme(panel.grid.major = element_blank(),
          legend.position = "none",
          axis.text.x=element_text(angle = 60, hjust = 1, size=input$cex_label_gene),
          axis.text.y=element_text(size=input$cex_label_category))
  if('fold.change' %in% names(data)){
    max.scale <- max(abs(res$fold.change))
    p <- p + scale_fill_distiller(type = "div",palette = "RdBu",
                                  direction = -1, limits = c(-max.scale, max.scale)) +
      theme(legend.position = "right")
  }
  return(p)
}