# Plot enrichment results (see run_gsea.R and run_ora.R)

plot_results<- function(resObject, geneList, ds.noext, output.dir, db.name, methodType){
  
  if(is.null(resObject))
    return()
  
  #trim descriptions to 80 characters
  resObject@result <- dplyr::mutate(resObject@result, Description = str_trunc(Description, 80))
  
  p <- enrichplot::dotplot(resObject, 
               showCategory = 20,
               label_format=50) 
  dotgr.fn <- paste(ds.noext, db.name, methodType, "dot_gr.pdf", sep = "_")
  ggplot2::ggsave(p, file = file.path(output.dir, methodType,"plots",dotgr.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  p <- enrichplot::dotplot(resObject, 
               showCategory = 20, 
               x = "count",
               label_format=50) 
  dotcnt.fn <- paste(ds.noext, db.name, methodType,"dot_cnt.pdf", sep = "_")
  ggplot2::ggsave(p, file = file.path(output.dir, methodType,"plots",dotcnt.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  # p <- enrichplot::cnetplot(resObject, foldChange=geneList,
  #               categorySize="geneNum", 
  #               cex_label_category = 0.8, 
  #               cex_label_gene = 1.0)
  # cnet1.fn <- paste(ds.noext, db.name, methodType,"cnet1.pdf", sep = "_")
  # ggplot2::ggsave(p, file = file.path(output.dir, methodType,"plots",cnet1.fn), 
  #        width = 2400, height = 2400, units = "px", device='pdf')
  # 
  # p <- enrichplot::cnetplot(resObject, foldChange=geneList, 
  #               circular = TRUE, colorEdge = TRUE, 
  #               cex_label_category = 0.8, 
  #               cex_label_gene = 1.0) 
  # cnet2.fn <- paste(ds.noext, db.name, methodType,"cnet2.pdf", sep = "_")
  # ggplot2::ggsave(p, file = file.path(output.dir, methodType,"plots",cnet2.fn), 
  #        width = 2400, height = 2400, units = "px", device='pdf')
  # 
  # data.emap <- pairwise_termsim(resObject)
  # p <- enrichplot::emapplot(data.emap, 
  #               showCategory = 20,
  #               cex_label_category=0.7,
  #               layout="nicely") #alt layouts: "kk","sugiyama","nicely","fr", "gem","lgl","mds"
  # emap.fn <- paste(ds.noext, db.name, methodType,"emap.pdf", sep = "_")
  # ggplot2::ggsave(p, file = file.path(output.dir, methodType,"plots",emap.fn), 
  #        width = 2400, height = 2400, units = "px", device='pdf')
  
  
  # # GSEA ONLY
  # if(methodType == "gsea"){
  #   p <- heatplot(resObject, foldChange=geneList, 
  #                 showCategory=5,
  #                 label_format=50) + coord_fixed(ratio=2)
  #   heat.fn <- paste(ds.noext, db.name, methodType,"heatmap.pdf", sep = "_")
  #   ggplot2::ggsave(p, file = file.path(output.dir, methodType,"plots",heat.fn), 
  #          width = 2400, height = 1200, units = "px", device='pdf')
  #   
  #   data.vol <- as.data.frame(resObject)
  #   p<-EnhancedVolcano(data.vol,
  #                      lab = data.vol$Description,
  #                      selectLab = head(data.vol$Description,3),
  #                      drawConnectors = TRUE,
  #                      widthConnectors = 0.2,
  #                      x = 'NES',
  #                      y = 'pvalue',
  #                      pCutoff = 0.05,
  #                      FCcutoff = 1,
  #                      legendLabels=c('NS','NES','p-value',
  #                                     'p-value & NES'),
  #                      xlab = 'Normalized Enrichment Score',
  #                      pointSize = 2.0,
  #                      labSize = 5.0)
  #   vol.fn <- paste(ds.noext, db.name, methodType,"volcano.pdf", sep = "_")
  #   ggplot2::ggsave(p, file = file.path(output.dir, methodType,"plots",vol.fn), 
  #          width = 2400, height = 2400, units = "px", device='pdf')
  #   
  #   for (i in 1:5){
  #     nes.fn <- paste(ds.noext, db.name, methodType,i,"nes.pdf", sep = "_")
  #     pdf(file.path(output.dir, methodType,"plots",nes.fn))
  #     p <- enrichplot::gseaplot2(resObject, 
  #                    geneSetID = i, 
  #                    title = resObject$Description[i]) 
  #     print(p)
  #     dev.off()
  #   }
  # }
  # 
  # # ORA only
  # if(methodType == "ora"){
  #   p <- heatplot(resObject, #can't figure out foldChange param
  #                 showCategory=5,
  #                 label_format=50) + coord_fixed(ratio=2)
  #   heat.fn <- paste(ds.noext, db.name, methodType,"heatmap.pdf", sep = "_")
  #   ggplot2::ggsave(p, file = file.path(output.dir, methodType,"plots",heat.fn), 
  #          width = 2400, height = 1200, units = "px", device='pdf')
  #   
  #   p <- upsetplot(resObject, n=10)
  #   upset.fn <- paste(ds.noext, db.name, methodType,"upset.pdf", sep = "_")
  #   ggplot2::ggsave(p, file = file.path(output.dir, methodType,"plots",upset.fn), 
  #          width = 2500, height = 1600, units = "px", device='pdf')
  #   
  #   data.vol <- as.data.frame(resObject)
  #   p<-EnhancedVolcano(data.vol,
  #                      lab = data.vol$Description,
  #                      selectLab = head(data.vol$Description,3),
  #                      drawConnectors = TRUE,
  #                      widthConnectors = 0.2,
  #                      x = 'Count',
  #                      y = 'pvalue',
  #                      pCutoff = 0.05,
  #                      FCcutoff = 5,
  #                      legendLabels=c('NS','Count','p-value',
  #                                     'p-value & Count'),
  #                      xlab = 'Number of Genes',
  #                      pointSize = 2.0,
  #                      labSize = 5.0)
  #   vol.fn <- paste(ds.noext, db.name, methodType,"volcano.pdf", sep = "_")
  #   ggplot2::ggsave(p, file = file.path(output.dir, methodType,"plots",vol.fn), 
  #          width = 2400, height = 2400, units = "px", device='pdf')
  # }
}