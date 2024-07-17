# Check R version
message("Checking R version")
x = readLines("http://cran.r-project.org/sources.html")
# the version number is in the next line of 'The latest release'
relver = gsub("(.*R-|\\.\\d\\.tar\\.gz.*)", "", x[grep("latest release", x) + 1L])
myver = substr(getRversion(),1,3)
# new version available?
# message("Installed major version: ", myver)
# message("Latest major version: ", relver)
vcheck <- compareVersion(relver, myver)
if(vcheck == 1){
  message(R.version$version.string)
  stop(sprintf("Please install R version %s or higher",relver))
}
# Check initial dependencies
message("Checking dependencies")
initial.libs <- c("shinydashboard", "shiny", "shinyjs", "shinyBS", "DT","BiocManager")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(initial.libs, try.bioconductor=TRUE, character.only = TRUE)
status <- sapply(initial.libs,require,character.only = TRUE)
  if(all(status)){
    print("All initial libraries successfully installed and loaded.")
  } else {   
     print(paste("ERROR: One or more libraries failed to install correctly.",
           status))
  }
# Additional dependencies loaded during analysis
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
  "rWikiPathways",
  "EnhancedVolcano")

# Install from cran 
p_load(p.load.libs, try.bioconductor=TRUE, character.only = TRUE)
status <- sapply(p.load.libs,require,character.only = TRUE)
if(all(status)){
  print("Additional libraries successfully installed and loaded.")
} else {   
  print(paste("ERROR: One or more libraries failed to install correctly.",
              status))
}
# Check which bioc packages are not installed
installed_packages <- installed.packages()[, "Package"]
packages_to_install <- bioc.load.libs[!(bioc.load.libs %in% installed_packages)]

# Install missing packages
if (length(packages_to_install) > 0) {
  BiocManager::install(packages_to_install, update = FALSE, force = TRUE)
} else {
  message("All specified packages are already installed.")
}


# Launch shiny app
shiny::runApp('shiny_result/app', launch.browser= TRUE)
