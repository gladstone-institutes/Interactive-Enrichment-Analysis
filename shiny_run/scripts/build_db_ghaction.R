# Process provided GMT files and build an RData collection

library(rWikiPathways)
library(dplyr)
library(magrittr)
library(jsonlite)
library(httr)

#GH Action version generates gmt.list and db.name at runtime
orgs <- c("Hs","Mm")
for(o in orgs){
  Sys.setenv(TZ="GMT")
  this.date <- format(Sys.Date(), "%Y%m%d")
  db.name <- paste(o,this.date, sep="_")
  gmt.list <- list.files("shiny_run/databases/gmts", pattern = paste0(".*_",tolower(o),"_.*.gmt"))
  
  #getFigureInfo for PFOCR figure filenames i/o aliases for jpg links
  r.info <- GET("https://pfocr.wikipathways.org/json/getFigureInfo.json")
  if(http_error(r.info)){
    stop("Failed to fetch PFOCR getFigureInfo.")
  }
  json_content <- content(r.info, as = "text", encoding = "UTF-8")
  json_list <- fromJSON(json_content, flatten = TRUE)
  df.pfocr <- as.data.frame(json_list) %>%
    dplyr::mutate(term = paste(figureInfo.pmcid,figureInfo.number,sep = "__")) %>%
    dplyr::rename(figid = figureInfo.figid) %>%
    dplyr::distinct(figid, term)
  
  #read gmt files and format gmt df objects
  for(g in gmt.list) {
    g.name <- gsub(".gmt","", g)
    g.fn <- file.path("shiny_run","databases","gmts",g)
    g.tg <- rWikiPathways::readGMT(g.fn)
    g.df <- g.tg
    # special handling
    if(grepl("%GOBP%", g.tg$term[1])){ #GO GMT
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
        dplyr::mutate(term = gsub(".jpg", "", term)) %>% #if it's there
        dplyr::left_join(df.pfocr, by=c("term"),relationship = "many-to-many") %>%
        dplyr::mutate(term = figid) %>%
        dplyr::select(term, name, gene)
    }
    
    assign(g.name, g.df)
  }
  
  db.name.list <- gsub(".gmt","", gmt.list)
  dir.create(file.path("shiny_run","databases","temp"))
  db.fn <- file.path("shiny_run","databases","temp",paste(db.name,"RData", sep = "."))
  save(list=db.name.list, file=db.fn)
}
