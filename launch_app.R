# Check R version
min.major = 4
min.minor = 2.0
R.version$version.string
if (R.version$major < min.major)
  stop(sprintf("ERROR: Please install R version %i.%f or higher",min.major, min.minor))
else 
  if (R.version$minor < min.minor)
    stop(sprintf("ERROR: Please install R version %i.%f or higher",min.major, min.minor))
    
# Check initial dependencies
initial.libs <- c("shinydashboard", "shiny", "shinyjs", "shinyBS", "DT","BiocManager")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(initial.libs, update = TRUE, try.bioconductor=TRUE, character.only = TRUE)
status <- sapply(initial.libs,require,character.only = TRUE)
  if(all(status)){
    print("All initial libraries successfully installed and loaded.")
  } else{   
     print(paste("ERROR: One or more libraries failed to install correctly.",
           status))
  }
# Launch shiny app
shiny::runGitHub('Interactive-Enrichment-Analysis', 'gladstone-institutes',subdir = 'shiny_run/app')
