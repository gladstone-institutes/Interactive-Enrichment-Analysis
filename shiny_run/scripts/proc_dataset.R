# Process dataset for enrichment analysis (see shiny_run)

proc_dataset<-function(ds.name, org.db.name, fromType, 
                       ora.fc, ora.pv, run.ora, run.gsea,
                       output.name="run"){
  
  # Objects from strings
  ds.fn <- file.path("../datasets",ds.name)
  dataset <- read.table(ds.fn, sep = ",", header = T, stringsAsFactors = F)
  names(dataset) <- tolower(names(dataset)) #force lowercase
  dataset$gene <- as.character(dataset$gene) #force characters
  ds.noext <- strsplit(ds.name,"\\.")[[1]][1]
  output.dir <- file.path("../",output.name, ds.noext)
  
  if(run.ora){ # Prepare subset of genes
    dir.create(file.path(output.dir,"ora","plots"), recursive = T, showWarnings = F)    
    set.genes <- NULL
    ds.names <- names(dataset)
    # TODO
    if(!'fold.change' %in% ds.names){
      if(!'p.value' %in% ds.names){ # plain list
        set.genes <- dataset%>%
          dplyr::mutate(ora.set = 1)
      } else { #p.value only
        set.genes <- dataset %>%
          rowwise() %>%
          dplyr::mutate(ora.set = ifelse(p.value < ora.pv, 1, 0)) %>%
          as.data.frame()
      }
    } else {
      if(!'p.value' %in% ds.names){
        return('Found "fold.change" but no "p.value." Please reformat your CSV to have a "p.value" column. Skipped!')
      } else { # fold.change and p.value 
        set.genes <- dataset %>%
          rowwise() %>%
          dplyr::mutate(ora.set = ifelse(p.value < ora.pv && abs(fold.change) > ora.fc , 1, 0)) %>%
          as.data.frame()
      }
    }
    
    # Identifier mapping
    set.genes.entrez <- map_ids(set.genes, org.db.name, fromType)
    
    # Record unmapped rows
    set.genes.unmapped <- set.genes %>%
      dplyr::filter(!gene %in% set.genes.entrez[[fromType]]) %>%
      {if('p.value' %in% ds.names) dplyr::arrange(.,p.value) else .} %>%
      dplyr::arrange(desc(ora.set))
    save_genes_params(set.genes.unmapped, ds.noext, "ora", output.dir, T, ora.fc, ora.pv)
    
    # Resolve duplicates (keep ENTREZID in ora.set and with smallest p.value, if available)
    set.genes.entrez.dedup <- set.genes.entrez %>%
      {if('p.value' %in% ds.names) dplyr::arrange(.,p.value) else .} %>%
      dplyr::arrange(desc(ora.set)) %>%
      dplyr::distinct(ENTREZID, .keep_all = T)
    save_genes_params(set.genes.entrez.dedup, ds.noext, "ora", output.dir, F, ora.fc, ora.pv)
  } 
  
  if (run.gsea){ # Rank list of genes
    dir.create(file.path(output.dir,"gsea","plots"), recursive = T, showWarnings = F)
    ranked.genes <- NULL
    ds.names <- names(dataset)
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
    ranked.genes.entrez <- map_ids(ranked.genes, org.db.name, fromType)
    
    # Record unmapped rows
    ranked.genes.unmapped <- ranked.genes %>%
      dplyr::filter(!gene %in% ranked.genes.entrez[[fromType]]) %>%
      dplyr::arrange(desc(rank))
    save_genes_params(ranked.genes.unmapped, ds.noext, "gsea", output.dir, T)
    
    # Resolve duplicates (keep ENTREZID with largest abs(rank))
    ranked.genes.entrez.dedup <- ranked.genes.entrez %>%
      dplyr::mutate(absrank = abs(rank)) %>%
      dplyr::arrange(desc(absrank)) %>%
      dplyr::distinct(ENTREZID, .keep_all = T) %>%
      dplyr::select(-absrank)  %>%
      dplyr::arrange(desc(rank))
    save_genes_params(ranked.genes.entrez.dedup, ds.noext, "gsea", output.dir, F)
  }
}

save_genes_params <- function(data, ds.noext, method.dir, output.dir, excluded=F, fc=1, pv=1e-05){
  #genes
  suffix <- "input"
  if (excluded)
    suffix <- "excluded"
  this.fn <- paste0(ds.noext, "__", method.dir, "_", suffix,".rds")
  saveRDS(data, file.path(output.dir,method.dir,this.fn))  
}

map_ids <- function(input,org.db.name, fromType){
  toType.list <- c("ENTREZID","SYMBOL")
  if (fromType == "ENTREZID") 
    toType.list <- c("SYMBOL")
  
  #perform ID mapping
  output <- bitr(input$gene, fromType = fromType,
                 toType = toType.list, 
                 OrgDb = eval(parse(text=org.db.name)))
  
  join_cols = c("gene")
  names(join_cols) <- fromType
  output <- dplyr::left_join(output, input, by=join_cols)  
  return(output)
}
