options(shiny.maxRequestSize=10*1024^2) #10MB file upload size limit

#Defaults and Parameters
minGSSize.default <- 3 # min database geneset size (see clusterProfiler)
maxGSSize.default <- 500 # max database geneset size (see clusterProfiler)
fc.default <- 1.0 # absolute value threshold for fold.change in ORA
pv.default <- 0.05 # threshold for p.value in ORA
analyzed.columns <- c('gene','fold.change','p.value', 'rank')

# editable list of supported orgs
supported.orgs <- list(human = "org.Hs.eg.db", 
                       mouse = "org.Mm.eg.db")

# editable list of supported idTypes (see keytypes per orgDb)
supported.idTypes <- c("SYMBOL",  
                       "ENSEMBL", 
                       "ENTREZID", 
                       "UNIPROT")

# editable list of named patterns for idTypes to support auto-detection,
# in order of specificity, more to less
supported.idTypes.patterns <- c(
  ENTREZID = "^\\d+$",
  ENSEMBL = "^ENS[A-Z]*[FPTG]\\d{11}$",
  UNIPROT = "^([A-N,R-Z][0-9][A-Z][A-Z, 0-9][A-Z, 0-9][0-9])|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\\.\\d+)?|([A-N,R-Z][0-9][A-Z][A-Z, 0-9][A-Z, 0-9][0-9][A-Z][A-Z, 0-9][A-Z, 0-9][0-9])$",
  SYMBOL = "[A-Za-z0-9]+"
)

# list of R libraries to load for analysis and visualization
p.load.libs <- c(
  "fs",
  "writexl",
  "stringr",
  "plyr",
  "tidyr",
  "dplyr",
  "magrittr",
  "rstudioapi",
  "ggplot2",
  "ggnewscale")
bioc.load.libs <- c(
  "AnnotationDbi",
  "DOSE",
  "enrichplot",
  "clusterProfiler",
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
# gmt.list <- list.files("../databases/gmts", ".gmt")
# if(length(gmt.list) > 1)
#   gmt.list <- c(gmt.list,"SELECT ALL")
# gmt.list <- c(gmt.list,"ADD NEW GMT FILES")
rdata.list <- list.files("../databases", ".RData")
rdata.list <- c(rdata.list,"BUILD NEW DATABASE")

# Check datasets
ds.list <- list.files("../datasets", "csv")
if(length(ds.list) > 1)
  ds.list <- c(ds.list,"SELECT ALL")
ds.list <- c(ds.list,"ADD NEW DATASETS")
