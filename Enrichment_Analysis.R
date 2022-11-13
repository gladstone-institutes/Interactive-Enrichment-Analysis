# main script to set configs and execute enrichment

library(fs)
source("scripts/load_libs.R")
source("scripts/proc_dataset.R")
source("scripts/run_gsea.R")
# source("run_ora.R")
source("scripts/setup.R")

#Defaults and Parameters
minGSSize <- 3 # see clusterProfiler
maxGSSize <- 500 # see clusterProfiler
ora.fc.default <- 1.0 # absolute value threshold for fold.change in ORA
ora.pv.default <- 0.05 # threshold for p.value in ORA
score.calculated <- TRUE 
output.name.default <- format(Sys.time(), "%Y%m%d_%H%M%S") # timestamp by default
supported.orgs <- list(human = "org.Hs.eg.db", 
                       mouse = "org.Mm.eg.db",
                       rat = "org.Rn.eg.db")

#Initialization
ora.fc <- NULL
ora.pv <- NULL
run.gsea <-'n'
run.ora <- 'n'

#Check and load
check_databases()
check_datasets()
p_load(org.db.name, update = TRUE, character.only = TRUE)

for (ds.name in ds.list){
  sprintf("Processing %s ", ds.name)
  ds.noext <- strsplit(ds.name,"\\.")[[1]][1]
  output.dir <- file.path(output.dir.root, ds.noext)
  geneList <- proc_dataset(ds.name, org.db.name, output.dir)
  for (db.name in db.list){
    if(run.gsea){
      sprintf("Running GSEA on %s using %s", ds.name, db.name)
      run_gsea(geneList, db.name, minGSSize, maxGSSize, 
               org.db.name, score.calculated, output.dir)
    }
    if(run.ora){
      sprintf("Running ORA on %s using %s", ds.name, db.name)
      #TODO
    }
  }
}
