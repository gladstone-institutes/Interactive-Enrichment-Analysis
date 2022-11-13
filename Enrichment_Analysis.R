# main script to set configs and execute enrichment

library(fs)
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
run.gsea <- 'n'
run.ora <- 'n'

#Check and load
check_database()
check_datasets()

for (db.name in db.list){
  for (ds.name in ds.list){
    if(run.gsea){
      sprintf("Running GSEA on %s using %s", ds.name, db.name)
      run_gsea(ds.name, db.name, minGSSize, maxGSSize, 
               org.db.name = org.db.name,
               score.calculated = score.calculated, 
               output.dir = output.dir)
    }
    if(run.ora){
      sprintf("Running ORA on %s using %s", ds.name, db.name)
      #TODO
    }
  }
}
