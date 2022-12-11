# Process provided GMT files and build an RData collection

library(rWikiPathways)
library(dplyr)
library(magrittr)

build_db <- function(gmt.list, db.name){
  
  #read gmt files and format gmt df objects
  for(g in gmt.list) {
    g.name <- gsub(".gmt","", g)
    g.fn <- file.path("..","databases","gmts",g)
    g.tg <- rWikiPathways::readGMT(g.fn)
    g.df <- g.tg
    # special handling
    if(startsWith(g.tg$term[1], "GOBP_")){ #GO GMT
      g.df <- g.tg %>%
        dplyr::mutate(term = gsub("GOBP_","",term)) %>%
        dplyr::mutate(term = gsub("_"," ",term)) %>%
        dplyr::mutate(term = tolower(term)) %>%
        dplyr::mutate(name = term)
    } else if (grepl("%WikiPathways_", g.tg$term[1])) { #WP GMT
      g.df <- rWikiPathways::readPathwayGMT(g.fn) %>%
        dplyr::rename(term = wpid) %>%
        dplyr::select(c(name, term, gene))

    } else if (grepl("^PMC\\d+__.*$", g.tg$term[1])) { # PFOCR
      g.df <- rWikiPathways::readGMTnames(g.fn) %>%
        dplyr::right_join(g.tg, by=c("term")) %>%
        dplyr::mutate(term = gsub(".jpg", "", term)) #if it's there
    }
    
    assign(g.name, g.df)
  }
  
  db.name.list <- gsub(".gmt","", gmt.list)
  db.fn <- file.path("..","databases",paste(db.name,"RData", sep = "."))
  save(list=db.name.list, file=db.fn)
  
  return()
}