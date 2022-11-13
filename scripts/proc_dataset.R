# Process dataset for enrichment analysis (see Enrichment_Analysis.R)

proc_dataset<-function(ds.name, db.name,output.dir="temp"){
  
  # Objects from strings
  ds.fn <- file.path("datasets",ds.name)
  dataset <- read.table(ds.fn, sep = ",", header = T, stringsAsFactors = F)
  ds.noext <- strsplit(ds.name,"\\.")[[1]][1]
  
  if(run.ora){ # Prepare split
    dir.create(file.path(output.dir,"ora","plots"), recursive = T, showWarnings = F)
    # TODO
    if(!'fold.change' %in% ds.names){
      if(!'p.value' %in% ds.names){
        if(!'rank' %in% ds.names){
        }
      }
    }
    
    p <- plot_volcano(ranked.genes.unmapped, ora.fc, ora.pv)
    vol1.fn <- paste(ds.noext,"volcano_excluded_unmapped.pdf", sep = "_")
    ggsave(p, file = file.path(output.dir,"ora","plots",vol1.fn), 
           width = 2400, height = 2400, units = "px", device='pdf')
    
    p <- plot_volcano(ranked.genes.entrez, ora.fc, ora.pv)
    vol2.fn <- paste(ds.noext,"volcano_all_genes.pdf", sep = "_")
    ggsave(p, file = file.path(output.dir,"ora","plots",vol2.fn), 
           width = 2400, height = 2400, units = "px", device='pdf')
    
    return(dataset) 
    
  } else if (run.gsea){ # Rank list of genes
    dir.create(file.path(output.dir,"gsea","plots"), recursive = T, showWarnings = F)
    ranked.genes <- NULL
    ds.names <- tolower(names(dataset))
    if(!'fold.change' %in% ds.names){
      if(!'p.value' %in% ds.names){
        if(!'rank' %in% ds.names){
          return('Could not find "rank" or "p.value" columns in this dataset. Skipped!')
        } else { #by user-provided rank
          ranked.genes <- dataset       
        }
      } else { #by user-provided p.value
        ranked.genes <- dataset %>%
          dplyr::mutate(rank = (-log10(as.numeric(as.character(p.value))))) 
      }
    } else {
      if(!'p.value' %in% ds.names){
        return('Found "fold.change" but no "p.value." Please reformat your CSV to have a "p.value" column. Skipped!')
      } else { #by user-provided fold.change and p.value
        ranked.genes <- dataset %>%
          dplyr::mutate(rank = (sign(as.numeric(fold.change)) * -log10(as.numeric(as.character(p.value))))) 
      }
    }
    
    # Identifier mapping
    ranked.genes.entrez <- bitr(ranked.genes$gene, fromType = "SYMBOL",
                                toType = c("ENTREZID"),
                                OrgDb = eval(parse(text=org.db.name)))
    ranked.genes.entrez <- dplyr::left_join(ranked.genes.entrez, ranked.genes, by = c("SYMBOL" = "gene"))
    ranked.genes.unmapped <- ranked.genes %>%
      dplyr::filter(!gene %in% ranked.genes.entrez$SYMBOL) %>%
      dplyr::arrange(desc(rank))
    #save and plot
    unmapped.fn <- paste0(ds.noext, "__excluded_unmapped.xlsx")
    write_xlsx(ranked.genes.unmapped, file.path(output.dir,"gsea",unmapped.fn))
    ranked.genes.unmapped$SYMBOL <- ranked.genes.unmapped$gene
    p <- plot_volcano(ranked.genes.unmapped)
    vol1.fn <- paste(ds.noext,"volcano_excluded_unmapped.pdf", sep = "_")
    ggsave(p, file = file.path(output.dir,"gsea","plots",vol1.fn), 
           width = 2400, height = 2400, units = "px", device='pdf')
    
    # Resolve duplicates (keep ENTREZID with largest abs(rank))
    ranked.genes.entrez.dedup <- ranked.genes.entrez %>%
      dplyr::mutate(absrank = abs(rank)) %>%
      dplyr::arrange(desc(absrank)) %>%
      dplyr::distinct(ENTREZID, .keep_all = T) %>%
      dplyr::select(-absrank)  %>%
      dplyr::arrange(desc(rank))
    #save and plot
    ranked.fn <- paste0(ds.noext, "__ranked_genes.xlsx")
    write_xlsx(ranked.genes.entrez.dedup, file.path(output.dir,"gsea",ranked.fn))
    p <- plot_volcano(ranked.genes.entrez.dedup)
    vol2.fn <- paste(ds.noext,"volcano_ranked_genes.pdf", sep = "_")
    ggsave(p, file = file.path(output.dir,"gsea","plots",vol2.fn), 
           width = 2400, height = 2400, units = "px", device='pdf')
    
    # Sorted named list for clusterProfiler, a.k.a. geneList
    ranked.genes.entrez.nl<-ranked.genes.entrez.dedup$rank
    names(ranked.genes.entrez.nl)<-ranked.genes.entrez.dedup$ENTREZID
    ranked.genes.entrez.nl <- sort(ranked.genes.entrez.nl, decreasing = T)
    
    return(ranked.genes.entrez.nl)
  }
}

plot_volcano<- function(data, fc=1, pv=1e-05){
  p<-EnhancedVolcano(data,
                     lab = data$SYMBOL,
                     selectLab = c(head(data$SYMBOL),tail(data$SYMBOL)),
                     x = 'fold.change',
                     y = 'p.value',
                     pCutoff = pv,
                     FCcutoff = fc,
                     legendLabels=c('NS','FC','p-value',
                                    'p-value & FC'),
                     pointSize = 2.0,
                     labSize = 5.0)
  return(p)
}