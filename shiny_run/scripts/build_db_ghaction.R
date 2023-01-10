# Process provided GMT files and build an RData collection

library(rWikiPathways)
library(dplyr)
library(magrittr)

#GH Action version generates gmt.list and db.name at runtime
orgs <- c("Hs","Mm")
for(o in orgs){
  Sys.setenv(TZ="GMT")
  this.date <- format(Sys.Date(), "%Y%m%d")
  db.name <- paste(o,this.date, sep="_")
  gmt.list <- list.files("shiny_run/databases/gmts", pattern = paste0(".*_",tolower(o),"_.*.gmt"))
  
  #read gmt files and format gmt df objects
  for(g in gmt.list) {
    g.name <- gsub(".gmt","", g)
    g.fn <- file.path("shiny_run","databases","gmts",g)
    g.tg <- rWikiPathways::readGMT(g.fn)
    g.df <- g.tg
    # special handling
    if(startsWith(g.tg$term[1], "GOBP_")){ #GO GMT
      g.df <- rWikiPathways::readPathwayGMT(g.fn) %>%
        dplyr::rename(term = wpid) %>%
        dplyr::mutate(name = paste0(toupper(substr(name,1,1)),tolower(substr(name,2,nchar(name))))) %>%
        dplyr::select(c(name, term, gene))
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
  dir.create(file.path("shiny_run","databases","temp"))
  db.fn <- file.path("shiny_run","databases","temp",paste(db.name,"RData", sep = "."))
  save(list=db.name.list, file=db.fn)
}
