## Function to perform GSEA (see Enrichment_Analysis.R)

run_gsea<-function(geneList, db.name, minGSSize, maxGSSize, org.db.name,
                   output.dir="temp"){
  
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
  gl.fn <- paste(ds.noext, db.name,"gsea","geneList.rds", sep = "_")
  saveRDS(geneList, file.path(output.dir,"gsea",gl.fn))
  er.fn <- paste(ds.noext, db.name,"gsea","result.rds", sep = "_")
  saveRDS(gseaResult, file.path(output.dir,"gsea",er.fn))
  
  ## Save df as TSV and XLSX
  gseaResult.df <- as.data.frame(gseaResult)
  tsv.fn <- paste(ds.noext, db.name,"gsea.tsv", sep = "_")
  xlsx.fn <- paste(ds.noext, db.name,"gsea.xlsx", sep = "_")
  write.table(gseaResult.df,file.path(output.dir,"gsea",tsv.fn),
              row.names=FALSE,sep="\t",quote=FALSE)
  write_xlsx(gseaResult.df,file.path(output.dir,"gsea",xlsx.fn))
  
  ## Plot
  plot_results(gseaResult, geneList, db.name, "gsea")
}
