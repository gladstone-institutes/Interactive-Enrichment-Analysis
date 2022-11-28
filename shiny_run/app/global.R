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
  "stringr",
  "dplyr",
  "magrittr",
  "clusterProfiler",
  "DOSE",
  "ggplot2",
  "ggupset",
  "enrichplot",
  "EnhancedVolcano")

jscode <- "
shinyjs.togglebox = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
}
shinyjs.hidebox = function(boxid){
$('#' + boxid).hide();
}
shinyjs.showbox = function(boxid){
$('#' + boxid).show();
}
shinyjs.startnew = function() {
history.go(0)
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
