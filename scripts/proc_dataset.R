# Process dataset for enrichment analysis (see Enrichment_Analysis.R)

proc_dataset<-function(ds.name, db.name,output.dir="temp"){
  
  # Objects from strings
  ds.fn <- file.path("datasets",ds.name)
  dataset <- read.table(ds.fn, sep = ",", header = T, stringsAsFactors = F)
  ds.noext <- strsplit(ds.name,"\\.")[[1]][1]
  
  if(run.ora){ # Prepare subset of genes
    dir.create(file.path(output.dir,"ora","plots"), recursive = T, showWarnings = F)    
    set.genes <- NULL
    ds.names <- tolower(names(dataset))
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
    set.genes.entrez <- map_ids(set.genes, fromType)
    
    # Record unmapped rows
    set.genes.unmapped <- set.genes %>%
      dplyr::filter(!gene %in% set.genes.entrez[[fromType]]) %>%
      {if('p.value' %in% ds.names) dplyr::arrange(.,p.value) else .} %>%
      dplyr::arrange(desc(ora.set))
    save_genes_params(set.genes.unmapped,"ora", T, ora.fc, ora.pv)
    
    # Resolve duplicates (keep ENTREZID in ora.set and with smallest p.value, if available)
    set.genes.entrez.dedup <- set.genes.entrez %>%
      {if('p.value' %in% ds.names) dplyr::arrange(.,p.value) else .} %>%
      dplyr::arrange(desc(ora.set)) %>%
      dplyr::distinct(ENTREZID, .keep_all = T)
    save_genes_params(set.genes.entrez.dedup, "ora", F, ora.fc, ora.pv)
    
    ora.geneList<<-set.genes.entrez.dedup
  } 
  
  if (run.gsea){ # Rank list of genes
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
    ranked.genes.entrez <- map_ids(ranked.genes, fromType)
    
    # Record unmapped rows
    ranked.genes.unmapped <- ranked.genes %>%
      dplyr::filter(!gene %in% ranked.genes.entrez[[fromType]]) %>%
      dplyr::arrange(desc(rank))
    save_genes_params(ranked.genes.unmapped, "gsea", T)
    
    # Resolve duplicates (keep ENTREZID with largest abs(rank))
    ranked.genes.entrez.dedup <- ranked.genes.entrez %>%
      dplyr::mutate(absrank = abs(rank)) %>%
      dplyr::arrange(desc(absrank)) %>%
      dplyr::distinct(ENTREZID, .keep_all = T) %>%
      dplyr::select(-absrank)  %>%
      dplyr::arrange(desc(rank))
    save_genes_params(ranked.genes.entrez.dedup, "gsea", F)
    
    # Sorted named list for clusterProfiler, a.k.a. geneList
    ranked.genes.entrez.nl<-ranked.genes.entrez.dedup$rank
    names(ranked.genes.entrez.nl)<-ranked.genes.entrez.dedup$ENTREZID
    ranked.genes.entrez.nl <- sort(ranked.genes.entrez.nl, decreasing = T)
    
    gsea.geneList<<-ranked.genes.entrez.nl
  }
}

save_genes_params <- function(data, method.dir, excluded=F, fc=1, pv=1e-05){
  #genes
  suffix <- "input"
  if (excluded)
    suffix <- "excluded"
  this.fn <- paste0(ds.noext, "__", method.dir, "_", suffix,".rds")
  saveRDS(data, file.path(output.dir,method.dir,this.fn))  
  #params
  par.df<-data.frame(fromType = fromType, pv=pv,fc=fc)
  this.fn <- paste0(ds.noext, "__", method.dir, "_params.rds")
  saveRDS(par.df, file.path(output.dir,method.dir,this.fn))  

  # 
  # # volcano plot if both p.value and fold.change are present
  # if('p.value' %in% names(data) & 'fold.change' %in% names(data)) {
  #   if(!fromType %in% names(data)) # i.e., unmapped cases
  #     data[fromType] <- data$gene
  #   
  #   p<-EnhancedVolcano(data,
  #                      lab = data[[fromType]],
  #                      selectLab = c(head(data[[fromType]]),tail(data[[fromType]])),
  #                      x = 'fold.change',
  #                      y = 'p.value',
  #                      pCutoff = pv,
  #                      FCcutoff = fc,
  #                      legendLabels=c('NS','FC','p-value',
  #                                     'p-value & FC'),
  #                      pointSize = 2.0,
  #                      labSize = 5.0)
  #   vol.fn <- paste0(ds.noext,"__",method.dir,"_volcano_",suffix,".pdf")
  #   ggsave(p, file = file.path(output.dir,method.dir,"plots",vol.fn), 
  #          width = 2400, height = 2400, units = "px", device='pdf')
  # }
}

map_ids <- function(input, fromType){
  output <- bitr(input$gene, fromType = fromType,
       toType = c("ENTREZID"),
       OrgDb = eval(parse(text=org.db.name)))
  join_cols = c("gene")
  names(join_cols) <- fromType
  output <- dplyr::left_join(output, input, by=join_cols)  
  return(output)
}
