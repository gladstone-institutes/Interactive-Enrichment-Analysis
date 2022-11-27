# Load files and variables

library(plyr)

#Defaults and Parameters
minGSSize.default <- 3 # min database geneset size (see clusterProfiler)
maxGSSize.default <- 500 # max database geneset size (see clusterProfiler)
fc.default <- 1.0 # absolute value threshold for fold.change in ORA
pv.default <- 0.05 # threshold for p.value in ORA
analyzed.columns <- c('gene','fold.change','p.value', 'rank')

# editable list of supported orgs
supported.orgs <- list(human = "org.Hs.eg.db", 
                       mouse = "org.Mm.eg.db",
                       rat = "org.Rn.eg.db")

# editable list of supported idTypes (see keytypes per orgDb)
supported.idTypes <- c("SYMBOL",  
                       "ENSEMBL", 
                       "ENTREZID", 
                       "UNIPROT")

# list of R libraries to load for analysis and visualization
load.libs <- c(
  "fs",
  "writexl",
  "dplyr",
  "magrittr",
  "clusterProfiler",
  "DOSE",
  "ggplot2",
  "ggupset",
  "enrichplot",
  "EnhancedVolcano")

# to collapse boxes
jscode <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
}
"

# Check databases
rdata.list <- list.files("../databases", ".RData")
rdata.list <- c(rdata.list,"BUILD NEW DATABASE")


# Check datasets
ds.list <- list.files("../datasets", "csv")
if(length(ds.list) > 1)
  ds.list <- c(ds.list,"SELECT ALL")
ds.list <- c(ds.list,"ADD NEW DATASETS")

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