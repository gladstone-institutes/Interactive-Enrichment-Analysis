## Function to perform ORA (see Enrichment_Analysis.R)

run_ora<-function(geneList, db.name, minGSSize, maxGSSize, org.db.name,
                   output.dir="temp"){
  
  # Object from string
  database <- eval(parse(text=db.name))
  
  gene <- geneList %>%
    dplyr::filter(ora.set == 1) %>%
    dplyr::pull(ENTREZID)
  universe <- pull(geneList, ENTREZID)
  
  # Use entire genome if not a susbset (i.e., by size or by missing p.value) 
  if(!'p.value' %in% names(geneList) | !length(universe) > length(gene)){
    universe <- NULL
  }
  
  # Perform ORA
  enrichResult <- clusterProfiler::enricher(
    gene,
    universe = universe,
    TERM2GENE = database[,c("term","gene")],
    TERM2NAME = database[,c("term","name")],
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    # pAdjustMethod="holm", #default is "BH"
    pvalueCutoff = 1 #to limit results
    )
  enrichResult <- setReadable(enrichResult, eval(parse(text=org.db.name)), keyType = "ENTREZID")
  
  # Save objects
  gl.fn <- paste(ds.noext, db.name,"oraGeneList.rds", sep = "_")
  saveRDS(gene, file.path(output.dir,"ora",gl.fn))
  er.fn <- paste(ds.noext, db.name,"enrichResult.rds", sep = "_")
  saveRDS(enrichResult, file.path(output.dir,"ora",er.fn))
  
  ## Save df as TSV and XLSX
  enrichResult.df <- as.data.frame(enrichResult)
  tsv.fn <- paste(ds.noext, db.name,"ora.tsv", sep = "_")
  xlsx.fn <- paste(ds.noext, db.name,"ora.xlsx", sep = "_")
  write.table(enrichResult.df,file.path(output.dir,"ora",tsv.fn),
              row.names=FALSE,sep="\t",quote=FALSE)
  write_xlsx(enrichResult.df,file.path(output.dir,"ora",xlsx.fn))
  
  ## Plot
  plot_results(enrichResult, gene, db.name, "ora")

}