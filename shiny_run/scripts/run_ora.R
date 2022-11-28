## Function to perform ORA (see Enrichment_Analysis.R)

run_ora<-function(dataset.name, db.name, output.name="run"){
  
  # file.prefix and output dir
  file.prefix <- strsplit(dataset.name,"\\.")[[1]][1] #remove ext if there
  output.dir <- file.path("../",output.name, file.prefix)
  
  # Retrieve params
  par.fn <- paste0(file.prefix, "__ora_params.rds")
  params <- readRDS(file.path(output.dir, "ora",par.fn))
  org.db.name <- params$org.db.name
  minGSSize <- params$minGSSize
  maxGSSize <- params$maxGSSize
  
  # Object from string
  database <- eval(parse(text=db.name))
  
  # Retrieve geneList 
  gl.fn <- paste0(file.prefix, "__ora_input.rds")
  geneList <- readRDS(file.path(output.dir, "ora",gl.fn))
  
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
  if(!is.null(enrichResult))
    enrichResult <- setReadable(enrichResult, eval(parse(text=org.db.name)), keyType = "ENTREZID")
  
  # Save objects
  gl.fn <- paste(file.prefix, db.name,"ora","geneList.rds", sep = "_")
  saveRDS(gene, file.path(output.dir,"ora",gl.fn))
  er.fn <- paste(file.prefix, db.name,"ora","result.rds", sep = "_")
  saveRDS(enrichResult, file.path(output.dir,"ora",er.fn))
  
  ## Save df as TSV and XLSX
  enrichResult.df <- as.data.frame(enrichResult)
  tsv.fn <- paste(file.prefix, db.name,"ora.tsv", sep = "_")
  xlsx.fn <- paste(file.prefix, db.name,"ora.xlsx", sep = "_")
  write.table(enrichResult.df,file.path(output.dir,"ora",tsv.fn),
              row.names=FALSE,sep="\t",quote=FALSE)
  write_xlsx(enrichResult.df,file.path(output.dir,"ora",xlsx.fn))
  
  ## Plot
  # plot_results(enrichResult, gene, db.name, "ora")

}