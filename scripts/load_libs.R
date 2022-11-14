## Load all library dependencies

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
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
  print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}