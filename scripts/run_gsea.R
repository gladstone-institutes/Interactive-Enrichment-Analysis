## Function to perform GSEA (see main.R)

load.libs <- c(
  "dplyr",
  "magrittr",
  "clusterProfiler",
  "DOSE",
  "GO.db",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "org.Rn.eg.db",
  "ggplot2",
  "ggupset",
  "enrichplot",
  "writexl")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
  print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}

run_gsea<-function(ds.name, db.name, minGSSize, maxGSSize, org.db.name,
                   score.calculated = FALSE, output.dir="temp"){
  
  # Objects from strings
  database <- eval(parse(text=db.name))
  ds.fn <- file.path("datasets",ds.name)
  dataset <- read.table(ds.fn, sep = ",", header = T, stringsAsFactors = F)
  ds.noext <- strsplit(ds.name,"\\.")[[1]][1]
  
  # Rank list of genes
  ranked.genes <- NULL
  ds.names <- tolower(names(dataset))
  if(!'fold.change' %in% ds.names){
    if(!'p.value' %in% ds.names){
      if(!'rank' %in% ds.names){
        return('Could not find "rank" or "p.value" columns in this dataset. Skipped!')
      } else { #by user-provided rank
        ranked.genes <- dataset %>%
          dplyr::arrange(desc(rank)) %>% #largest values first
          dplyr::select(gene,rank)
      }
    } else { #by user-provided p.value
      ranked.genes <- dataset %>%
        dplyr::mutate(rank = (-log10(as.numeric(as.character(p.value))))) %>%
        dplyr::arrange(desc(rank)) %>% #largest -log10(pv) first (all positive; only top of list is of interest)
        dplyr::select(gene,rank)
    }
  } else {
    if(!'p.value' %in% ds.names){
      return('Found "fold.change" but no "p.value." Please reformat your CSV to have a "p.value" column. Skipped!')
    } else { #by user-provided fold.change and p.value
      ranked.genes <- dataset %>%
        dplyr::mutate(rank = (sign(as.numeric(fold.change)) * -log10(as.numeric(as.character(p.value))))) %>%
        dplyr::arrange(desc(rank)) %>% #largest positive values first (both ends of list are of interest)
        dplyr::select(gene,rank)
    }
  }
    
  # Identifier mapping
  ranked.genes.entrez <- bitr(ranked.genes$gene, fromType = "SYMBOL",
                                   toType = c("ENTREZID"),
                                   OrgDb = eval(parse(text=org.db.name)))
  
  ranked.genes.entrez <- dplyr::left_join(ranked.genes.entrez, ranked.genes, by = c("SYMBOL" = "gene"))
  ranked.genes.entrez.nl<-ranked.genes.entrez$rank
  names(ranked.genes.entrez.nl)<-ranked.genes.entrez$ENTREZID
  
  enrichment.result <- GSEA(
    ranked.genes.entrez.nl,
  #  pAdjustMethod="holm", #default is "BH"
    TERM2GENE = database[,c("term","gene")],
    TERM2NAME = database[,c("term","name")],
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = 1, #for results
    verbose=FALSE)
  enrichment.result <- setReadable(enrichment.result, eval(parse(text=org.db.name)), keyType = "ENTREZID")

  # Saving results
  dir.create(file.path(output.dir,"gsea"), showWarnings = F)
  dir.create(file.path(output.dir,"gsea","plots"), showWarnings = F)
  
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
  
  p <- cnetplot(enrichment.result, foldChange=ranked.genes.entrez.nl,
           categorySize="geneNum", 
           cex_label_category = 0.8, 
           cex_label_gene = 1.0)
  cnet1.fn <- paste(ds.noext, db.name,"cnet1.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",cnet1.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  p <- cnetplot(enrichment.result, foldChange=ranked.genes.entrez.nl, 
           circular = TRUE, colorEdge = TRUE, 
           cex_label_category = 0.8, 
           cex_label_gene = 1.0) 
  cnet2.fn <- paste(ds.noext, db.name,"cnet2.pdf", sep = "_")
  ggsave(p, file = file.path(output.dir,"gsea","plots",cnet2.fn), 
         width = 2400, height = 2400, units = "px", device='pdf')
  
  p <- heatplot(enrichment.result, foldChange=ranked.genes.entrez.nl, 
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
  enrichment.result.df <- as.data.frame(enrichment.result)
  tsv.fn <- paste(ds.noext, db.name,".tsv", sep = "_")
  xlsx.fn <- paste(ds.noext, db.name,".xlsx", sep = "_")
  write.table(enrichment.result.df,file.path(output.dir,"gsea",tsv.fn),
              row.names=FALSE,sep="\t",quote=FALSE)
  write_xlsx(enrichment.result.df,file.path(output.dir,"gsea",xlsx.fn))
  

  return()
}
