# Load files and variables

library(plyr)

#Defaults and Parameters
minGSSize <- 3 # min database geneset size (see clusterProfiler)
maxGSSize <- 500 # max database geneset size (see clusterProfiler)
ora.fc.default <- 1.0 # absolute value threshold for fold.change in ORA
ora.pv.default <- 0.05 # threshold for p.value in ORA
output.name <- format(Sys.time(), "%Y%m%d_%H%M%S") # timestamp 
supported.orgs <- list(human = "org.Hs.eg.db", # editable list of supported orgs
                       mouse = "org.Mm.eg.db",
                       rat = "org.Rn.eg.db")
supported.idTypes <- c("SYMBOL",  # editable list of supported idTypes (see keytypes per orgDb)
                       "ENSEMBL", 
                       "ENTREZID", 
                       "UNIPROT")

# Check databases
rdata.list <- list.files("../databases", ".RData")
rdata.list <- c(rdata.list,"Build a new database from GMTs")


# Check datasets
ds.list <- list.files("../datasets", "csv")
ds.list <- c(ds.list,"Add a dataset")

# # Dataframe of results in output dir
# output.dirs <- list.dirs("../output", full.names = F, recursive = T)
# output.df <- plyr::ldply(output.dirs, function(d){
#   d0 <- strsplit(d, "\\/")[[1]]
#   if(length(d0) == 4)
#     return(d0)
# }) 
# 
# if(nrow(output.df) < 1)
#   stop("No results found in output folder. Please run analysis first.")
# 
# names(output.df) <- c("run","dataset","method","sub")
# 
# 
# # Initialize lists
# run.list <- unique(output.df$run)
# ds.list <- unique(output.df[which(output.df$run==run.list[1]),'dataset'])
# method.list <- unique(output.df[which(output.df$run==run.list[1] & 
#                                         output.df$dataset==ds.list[1]),'method'])
# 
# db.list <- c("go_20210108", "wp_20210108", "pfocr_20210108")
# 


# #place the RDS file for the single-cell data in the data folder
# #update the name of the RDS file in the below line
# data <- readRDS(file = 'data/example_data.rds')
# #update the prefix for the filename of any data/plot downloaded from the app
# plot_download_prefix <- "example_app"