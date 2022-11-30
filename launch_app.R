# Verify initial dependencies
install.packages(c("shinydashboard", "shiny", "shinyjs", "shinyBS", "DT","BiocManager"))

# Launch shiny app
shiny::runGitHub('Interactive-Enrichment-Analysis', 'gladstone-institutes',subdir = 'shiny_run/app')