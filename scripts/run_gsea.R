## Function to perform GSEA (see Enrichment_Analysis.R)

run_gsea<-function(ds.proc, db.name, minGSSize, maxGSSize, org.db.name,
                   score.calculated = FALSE, output.dir="temp"){
  
  # Object from string
  database <- eval(parse(text=db.name))

  # Perform GSEA
  enrichment.result <- clusterProfiler::GSEA(
    ds.proc,
    TERM2GENE = database[,c("term","gene")],
    TERM2NAME = database[,c("term","name")],
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    # pAdjustMethod="holm", #default is "BH"
    pvalueCutoff = 1, #to limit results
    verbose=FALSE)
  enrichment.result <- setReadable(enrichment.result, eval(parse(text=org.db.name)), keyType = "ENTREZID")
  
  ## Plots
  p <- dotplot(enrichment.result, 
          showCategory = 20,
          label_format=50) 
  dotgr.fn <- paste(ds.noext, db.name,"dot_gr.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",dotgr.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  p <- dotplot(enrichment.result, 
          showCategory = 20, 
          x = "count",
          label_format=50) 
  dotcnt.fn <- paste(ds.noext, db.name,"dot_cnt.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",dotcnt.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  p <- cnetplot(enrichment.result, foldChange=ds.proc,
           categorySize="geneNum", 
           cex_label_category = 0.8, 
           cex_label_gene = 1.0)
  cnet1.fn <- paste(ds.noext, db.name,"cnet1.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",cnet1.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  p <- cnetplot(enrichment.result, foldChange=ds.proc, 
           circular = TRUE, colorEdge = TRUE, 
           cex_label_category = 0.8, 
           cex_label_gene = 1.0) 
  cnet2.fn <- paste(ds.noext, db.name,"cnet2.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",cnet2.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  p <- heatplot(enrichment.result, foldChange=ds.proc, 
           showCategory=10,
           label_format=50) + coord_fixed(ratio=2)
  heat.fn <- paste(ds.noext, db.name,"heatmap.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",heat.fn), 
         width = 2400, height = 1200, units = "px", device='pdf')
  
  enrichment.result.emap <- pairwise_termsim(enrichment.result)
  p <- emapplot(enrichment.result.emap, 
           showCategory = 20,
           cex_label_category=0.7,
           layout="nicely") #alt layouts: "kk","sugiyama","nicely","fr", "gem","lgl","mds"
  emap.fn <- paste(ds.noext, db.name,"emap.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",emap.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  p <- upsetplot(enrichment.result, n=10)
  upset.fn <- paste(ds.noext, db.name,"upset.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",upset.fn), 
         width = 2500, height = 1600, units = "px", device='pdf')
  
  enrichment.result.df <- as.data.frame(enrichment.result)
  p<-EnhancedVolcano(enrichment.result.df,
                     lab = enrichment.result.df$Description,
                     selectLab = head(enrichment.result.df$Description,3),
                     drawConnectors = TRUE,
                     widthConnectors = 0.2,
                     x = 'NES',
                     y = 'pvalue',
                     pCutoff = 1e-05,
                     FCcutoff = 1,
                     legendLabels=c('NS','NES','p-value',
                                    'p-value & NES'),
                     xlab = 'Normalized Enrichment Score',
                     pointSize = 2.0,
                     labSize = 5.0)
  vol.fn <- paste(ds.noext, db.name,"volcano.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",vol.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  for (i in 1:5){
    nes.fn <- paste(ds.noext, db.name,i,"nes.pdf", sep = "_")
    pdf(file.path(output.dir,"gsea","plots",nes.fn))
    p <- gseaplot2(enrichment.result, 
                  geneSetID = i, 
                  title = enrichment.result$Description[i]) 
    print(p)
    dev.off()
    # ggsave(pgrid, file = file.path(output.dir,"gsea","plots",nes.fn), 
    #        width = 2500, height = 2000, units = "px", device='pdf')
  }
  
  ## Write to TSV and XLSX
  tsv.fn <- paste(ds.noext, db.name,"gsea.tsv", sep = "_")
  xlsx.fn <- paste(ds.noext, db.name,"gsea.xlsx", sep = "_")
  write.table(enrichment.result.df,file.path(output.dir,"gsea",tsv.fn),
              row.names=FALSE,sep="\t",quote=FALSE)
  write_xlsx(enrichment.result.df,file.path(output.dir,"gsea",xlsx.fn))
  

  return()
}
