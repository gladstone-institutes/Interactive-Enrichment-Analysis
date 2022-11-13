# Setup functions for enrichment analysis (see main.R)

check_database <- function(update=FALSE){
  if (!update & exists("db.list")){
    print("These are the databases we will use in this run:")
    return(print.data.frame(as.data.frame(db.list)))
    check_db_colnames()
  }
  rdata.list <- list.files("database", ".RData")
  if(length(rdata.list) == 0){
    stop("No RData found in database folder. Please run prep_database() prior to run.")
    #TODO: check_gmts() to read and process gmt files in RData and then load_databases()
  } else if (length(rdata.list) == 1){
    rdata.fn <- file.path("database",rdata.list[1])
    load_database(rdata.fn)
  } else {
    print("Which database object should we load for this run:")
    print.data.frame(as.data.frame(rdata.list))
    rdata.pick <- readline("Enter row number: ")
    rdata.fn <- file.path("database",rdata.list[as.numeric(rdata.pick)])
    load_database(rdata.fn)
  }
}

load_database <- function(rdata.fn){
  db.list <<- load(rdata.fn, globalenv())
  print("These are the databases we will use in this run:")
  print.data.frame(as.data.frame(db.list))
  check_db_colnames()
  check_db_long_names()
}

check_db_colnames <- function(){
  for (db in db.list){
    db.colnames <- tolower(names(eval(parse(text=db))))
    for (cn in c("name","term","gene")){
      if(!cn %in% db.colnames){
        stop(sprintf('%s is missing a %s column',db, cn))
      }
    }
  }
}

check_db_long_names<- function(){
  for (db in db.list){
    db.df <- eval(parse(text=db)) %>%
      dplyr::mutate(name = stringr::str_trunc(name, 87))
    assign(db, db.df, envir = globalenv())
  }
}

check_datasets <- function(){
  ds.list <<- list.files("datasets", ".csv")
  if(length(ds.list) == 0){
    stop("No CSV found in datasets folder. Please format CSV files for analysis.")
  } 
  print("These are the datasets we will analyze in this run:")
  print.data.frame(as.data.frame(ds.list))
  
  print("Assuming all datasets are of the same format, let's examine the first one...")
  ds.fn <- file.path("datasets",ds.list[1])
  ds.one <- read.table(ds.fn, sep = ",", header = T, stringsAsFactors = F)
  
  print(head(ds.one))
  
  ds.names <- tolower(names(ds.one))
  #GENE CHECK
  if(!'gene' %in% ds.names){
    stop('Please reformat your CSV to have a "gene" column with gene names.')
  }
  #RANK AND THRESHOLD CHECKS
  if(!'fold.change' %in% ds.names){
    if(!'p.value' %in% ds.names){
      if(!'rank' %in% ds.names){
        run.ora <<- readline('Could not find "rank" or "p.value" columns. Proceed with ORA with the entire genome as background. (Y/n): ')
      } else {
        run.gsea <<- readline('Found "rank" column. Proceed with GSEA analysis? (Y/n): ')
      }
    } else {
      run.gsea <<- readline('Found "p.value" column. Proceed with GSEA analysis? (Y/n): ')
      run.ora <<- readline('Also perform ORA analysis? (Y/n): ')
      if (tolower(run.ora) %in% c("","y")){
        ora.pv <<- as.numeric(readline('  What threshold for p.value? (Hit Return for default of 0.05 or enter value): '))
      }
    }
  } else {
    if(!'p.value' %in% ds.names){
      stop('Found "fold.change" but no "p.value." Please reformat your CSV to have a "p.value" column.')
    } else {
      run.gsea <<- readline('Found "fold.change" and "p.value" columns. Proceed with GSEA analysis? (Y/n): ')
      run.ora <<- readline('Also perform ORA analysis? (Y/n): ')
      if (tolower(run.ora) %in% c("","y")){
        ora.fc <<- as.numeric(readline('  What threshold for absolute value of fold.change? (Hit Return for default of 1.0 or enter value): '))
        ora.pv <<- as.numeric(readline('  What threshold for p.value? (Hit Return for default of 0.05 or enter value): '))
      }
    }
  }
  
  #Process answers
  if (tolower(run.ora) %in% c("","y"))
    run.ora <<- TRUE
  else
    run.ora <<- FALSE
  if (tolower(run.gsea) %in% c("","y"))
    run.gsea <<- TRUE
  else
    run.gsea <<- FALSE
  if (!is.null(ora.fc))
    if (is.na(ora.fc)) #non-numeric
      ora.fc <<- ora.fc.default
  if (!is.null(ora.pv))
    if (is.na(ora.pv))#non-numeric
      ora.pv <<- ora.pv.default
  
  # Organism
  cat('\n')
  org.list <- names(supported.orgs)
  print("Which organism are we working with?")
  print.data.frame(as.data.frame(org.list))
  org.pick <- as.numeric(readline("Enter row number (or hit Return for human): "))
  if (is.na(org.pick)) #non-numeric response
    org.pick <- 1
  org.name <<- org.list[as.numeric(org.pick)]
  org.db.name <<- unlist(unname(supported.orgs[org.name]))
  
  # Run title
  cat('\n')
  output.name <<- readline('Custom name for this run? (Hit Return to use timestamp): ')
  if (output.name == "")
    output.name <<- output.name.default
  
  output.name <<- fs::path_sanitize(output.name, replacement = "")
  output.name <<- gsub("^\\W+", "", output.name) #clean up any lingering non-alnum at start
  
  safe.dir <- 'y'
  output.dir <<- file.path("output",output.name)
  if(dir.exists(output.dir))
    safe.dir <- readline('Output folder with the same name already exists. Overwrite files? (y/N): ')
  if(tolower(safe.dir) == 'y'){
    dir.create(output.dir, showWarnings = F)
    sprintf('Analysis results will be found in %s/', output.dir)
  } else {
    stop('Please try again. A unique name for the run is required (or overwrite confirmation).')
  }
}
