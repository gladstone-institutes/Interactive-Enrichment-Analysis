## Function to perform GSEA (see Enrichment_Analysis.R)

run_gsea<-function(geneList, db.name, minGSSize, maxGSSize, org.db.name,
                   score.calculated = FALSE, output.dir="temp"){
  
  # Object from string
  database <- eval(parse(text=db.name))

  # Perform GSEA
  gseaResult <- clusterProfiler::GSEA(
    geneList,
    TERM2GENE = database[,c("term","gene")],
    TERM2NAME = database[,c("term","name")],
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    # pAdjustMethod="holm", #default is "BH"
    pvalueCutoff = 1, #to limit results
    verbose=FALSE)
  gseaResult <- setReadable(gseaResult, eval(parse(text=org.db.name)), keyType = "ENTREZID")

  # Save objects
  gl.fn <- paste(ds.noext, db.name,"geneList.rds", sep = "_")
  saveRDS(geneList, file.path(output.dir,"gsea",gl.fn))
  er.fn <- paste(ds.noext, db.name,"gseaResult.rds", sep = "_")
  saveRDS(gseaResult, file.path(output.dir,"gsea",er.fn))
  
  ## Save df as TSV and XLSX
  gseaResult.df <- as.data.frame(gseaResult)
  tsv.fn <- paste(ds.noext, db.name,"gsea.tsv", sep = "_")
  xlsx.fn <- paste(ds.noext, db.name,"gsea.xlsx", sep = "_")
  write.table(gseaResult.df,file.path(output.dir,"gsea",tsv.fn),
              row.names=FALSE,sep="\t",quote=FALSE)
  write_xlsx(gseaResult.df,file.path(output.dir,"gsea",xlsx.fn))
  
  ## Plot
  # TODO: generate dynamically in RShiny
  # only vars needed (gseaResult, geneList) are saved as RDS
  p <- dotplot(gseaResult, 
          showCategory = 20,
          label_format=50) 
  dotgr.fn <- paste(ds.noext, db.name,"dot_gr.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",dotgr.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  p <- dotplot(gseaResult, 
          showCategory = 20, 
          x = "count",
          label_format=50) 
  dotcnt.fn <- paste(ds.noext, db.name,"dot_cnt.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",dotcnt.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  p <- cnetplot(gseaResult, foldChange=geneList,
           categorySize="geneNum", 
           cex_label_category = 0.8, 
           cex_label_gene = 1.0)
  cnet1.fn <- paste(ds.noext, db.name,"cnet1.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",cnet1.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  p <- cnetplot(gseaResult, foldChange=geneList, 
           circular = TRUE, colorEdge = TRUE, 
           cex_label_category = 0.8, 
           cex_label_gene = 1.0) 
  cnet2.fn <- paste(ds.noext, db.name,"cnet2.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",cnet2.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  p <- heatplot(gseaResult, foldChange=geneList, 
           showCategory=10,
           label_format=50) + coord_fixed(ratio=2)
  heat.fn <- paste(ds.noext, db.name,"heatmap.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",heat.fn), 
         width = 2400, height = 1200, units = "px", device='pdf')
  
  data.emap <- pairwise_termsim(gseaResult)
  p <- emapplot(data.emap, 
           showCategory = 20,
           cex_label_category=0.7,
           layout="nicely") #alt layouts: "kk","sugiyama","nicely","fr", "gem","lgl","mds"
  emap.fn <- paste(ds.noext, db.name,"emap.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",emap.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  p <- upsetplot(gseaResult, n=10)
  upset.fn <- paste(ds.noext, db.name,"upset.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",upset.fn), 
         width = 2500, height = 1600, units = "px", device='pdf')
  
  data.vol <- as.data.frame(gseaResult)
  p<-EnhancedVolcano(data.vol,
                     lab = data.vol$Description,
                     selectLab = head(data.vol$Description,3),
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
    p <- gseaplot2(gseaResult, 
                  geneSetID = i, 
                  title = gseaResult$Description[i]) 
    print(p)
    dev.off()
  }
  
  return()
}
