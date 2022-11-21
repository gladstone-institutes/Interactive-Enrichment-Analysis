# Load files and variables
library(plyr)

#Dataframe of results in output dir
output.dirs <- list.dirs("../output", full.names = F, recursive = T)
output.df <- plyr::ldply(output.dirs, function(d){
  d0 <- strsplit(d, "\\/")[[1]]
  if(length(d0) == 4)
    return(d0)
}) 

if(nrow(output.df) < 1)
  stop("No results found in output folder. Please run analysis first.")

names(output.df) <- c("run","dataset","method","sub")

# Initialize lists
run.list <- unique(output.df$run)
ds.list <- unique(output.df[which(output.df$run==run.list[1]),'dataset'])
method.list <- unique(output.df[which(output.df$run==run.list[1] & 
                                        output.df$dataset==ds.list[1]),'method'])

# #place the RDS file for the single-cell data in the data folder
# #update the name of the RDS file in the below line
# data <- readRDS(file = 'data/example_data.rds')
# #update the prefix for the filename of any data/plot downloaded from the app
# plot_download_prefix <- "example_app"