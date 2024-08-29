# Check R version
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
initial.libs <- c("shinydashboard", "shiny", "shinyjs", "shinyBS", 
                  "rstudioapi", "DT","BiocManager", "devtools",
                  "httr")
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

# Install or launch shiny app
local_dir <- file.path(getwd(), "Interactive-Enrichment-Analysis")

# Check if already installed
if (!dir.exists(local_dir)) {
  message("Installation not found. Downloading and running the app...")
  shiny::runGitHub('Interactive-Enrichment-Analysis', 
                   'gladstone-institutes',
                   subdir = 'shiny_run/app', 
                   destdir = local_dir, 
                   launch.browser = TRUE)
} else {
  message("Installation found. Launching the app...")
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
  app_dir <- file.path(script_dir, 'shiny_run/app')
  shiny::runApp(app_dir, launch.browser = TRUE)
}
