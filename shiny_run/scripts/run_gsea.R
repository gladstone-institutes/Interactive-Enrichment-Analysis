## Function to perform GSEA (see shiny_run)

run_gsea<-function(dataset.name, db.name, output.name="run"){
  
  # file.prefix and output dir
  file.prefix <- strsplit(dataset.name,"\\.")[[1]][1] #remove ext if there
  output.dir <- file.path("../",output.name, file.prefix)
  
  # Retrieve params
  par.fn <- paste0(file.prefix, "__gsea_params.rds")
  params <- readRDS(file.path(output.dir, "gsea",par.fn))
  org.db.name <- params$org.db.name
  minGSSize <- params$minGSSize
  maxGSSize <- params$maxGSSize
  
  # Object from string
  database <- get(db.name)
  
  # geneList from file.prefix
  gl.fn <- paste0(file.prefix, "__gsea_input.rds")
  geneList <- readRDS(file.path(output.dir, "gsea",gl.fn))

  # Sorted named list for clusterProfiler, a.k.a. geneList
  ## pre-ranking
  # geneList <- geneList %>%
  #   mutate(rank = rank(rank,  ties.method = "random")) %>%
  #   arrange(desc(rank))
  ranked.genes.entrez.nl<-geneList$rank
  names(ranked.genes.entrez.nl)<-geneList$ENTREZID
  geneList <- sort(ranked.genes.entrez.nl, decreasing = T)
  
  # Perform GSEA
  gseaResult <- clusterProfiler::GSEA(
    geneList,
    TERM2GENE = database[,c("term","gene")],
    TERM2NAME = database[,c("term","name")],
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    # pAdjustMethod="holm", #default is "BH"
    pvalueCutoff = 0.1, #to limit results
    verbose=FALSE)
  if(!is.null(gseaResult))
    gseaResult <- setReadable(gseaResult, get(org.db.name), keyType = "ENTREZID")
  
  # Save objects
  gl.fn <- paste(file.prefix, db.name,"gsea","geneList.rds", sep = "_")
  saveRDS(geneList, file.path(output.dir,"gsea",gl.fn))
  er.fn <- paste(file.prefix, db.name,"gsea","result.rds", sep = "_")
  saveRDS(gseaResult, file.path(output.dir,"gsea",er.fn))
  
  ## Save df as TSV and XLSX
  gseaResult.df <- as.data.frame(gseaResult)
  tsv.fn <- paste(file.prefix, db.name,"gsea.tsv", sep = "_")
  xlsx.fn <- paste(file.prefix, db.name,"gsea.xlsx", sep = "_")
  write.table(gseaResult.df,file.path(output.dir,"gsea",tsv.fn),
              row.names=FALSE,sep="\t",quote=FALSE)
  write_xlsx(gseaResult.df,file.path(output.dir,"gsea",xlsx.fn))
  
  ## Plot
  if (nrow(gseaResult.df  ) > 0)
    plot_results(gseaResult, geneList, file.prefix, output.dir, db.name, "gsea")
}
