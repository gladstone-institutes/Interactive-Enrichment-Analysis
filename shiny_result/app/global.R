library(plyr)

# Dataframe of results in output dir
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
run.default <- tail(run.list, 1)
ds.list <- unique(output.df[which(output.df$run==run.default),'dataset'])
method.list <- rev(unique(output.df[which(output.df$run==run.default & 
                                        output.df$dataset==ds.list[1]),'method']))
par.fp <- file.path("../output",run.default, ds.list[1], method.list[1], 
                paste0(paste(ds.list[1], method.list[1], sep = "__"), "_params.rds"))
par.df <- readRDS(par.fp)
db.list <- strsplit(par.df$db.list, ",")[[1]]
# db.list <- c("go_20210108", "wp_20210108", "pfocr_20210108")
