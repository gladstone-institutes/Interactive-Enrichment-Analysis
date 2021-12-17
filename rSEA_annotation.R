
database_lists1<-load("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/databases.RData")#has wp, pfocr, go
database_lists1=c(database_lists1,load("/Users/mingyoungshin/Dropbox (Gladstone)/enrichment_databases/jensen_database.RData"))

min_set_size <- 3
max_set_size <- 500 

database_lists1 <- unname(unlist(sapply(database_lists1, grep, pattern="_list$", value = T, perl = T)))
for (db in database_lists1) {
  eval(call("<-", as.name(db),  Filter(Negate(is.null), lapply(get(db), function(x){
    if(length(x) < min_set_size | length(x) > max_set_size)
      NULL
    else
      x
  }))
  ))
}


require(rSEA)
#data_m=down;database=go_list;database_annotation=go_annotation;output_name="down_GO_test"

get_rsea<-function(data_m,database,database_annotation, output_name="temp"){
  #merged=merge(data_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
  #merged=merged[!is.na(merged$ENTREZID),]
  
  rsea_result<-SEA(as.numeric(as.character(data_m$p.value)), data_m$ENTREZID, pathlist = database)
  

  result=NULL

    for(j in 1:nrow(rsea_result)){
      annotation=database_annotation[which(database_annotation$set_id==rsea_result[j,2])[1],1]
      genes=database_annotation[which(database_annotation$set_id==rsea_result[j,2]),3]
      genes=data_m$SYMBOL[which(as.character(data_m$ENTREZID)%in%as.character(genes))]
      #genes=genes[which(genes%in%data_m$Gene)]
      all_genes=paste(genes,collapse="/")
    
      
      gene_hits=database_annotation[which(database_annotation$set_id==rsea_result[j,2]),3]
      gene_hits=data_m$SYMBOL[intersect(which(as.character(data_m$ENTREZID)%in%as.character(gene_hits))
                              ,which(data_m$q.value<0.05) )]
      
      gene_hits=paste(gene_hits,collapse="/")
      
      result=rbind(result,cbind(rsea_result[j,],annotation,all_genes,gene_hits))
      
    }
    
    result=result[order(result$Comp.adjP,decreasing = FALSE),]
    result=result[!is.na(result[,5]),]
    write.table(result,paste0(output_name,"_all_results.txt"),col.names=colnames(result),row.names=FALSE,sep="\t",quote=FALSE)
    
    
 
  
  
}



match_entrez<-function(data){
  gene_entrez <- clusterProfiler::bitr(data$SYMBOL, fromType = "SYMBOL",
                                       toType = c("ENTREZID"),
                                       OrgDb = org.Hs.eg.db, drop=FALSE)
  
  tb=table(gene_entrez$SYMBOL)
  print(max(tb[which(tb>1)]))
  
  if(max(tb[which(tb>1)])>1){
    
    remove_i=NULL
    for(n in names(tb[which(tb>1)])){
      all_i=which(n==gene_entrez$SYMBOL)
      remove_i=c(remove_i,all_i[-1])
    }
    
    print(dim(gene_entrez))
    gene_entrez=gene_entrez[-remove_i,]
  }

  
  
  print(dim(gene_entrez))
  print(dim(data))
  data=merge(data,gene_entrez,by="SYMBOL")
  print(dim(data))
  colnames(data)[ncol(data)]="ENTREZID"
  
  print(dim(data))
  return(data)
}

library(magrittr)
library(org.Hs.eg.db)

#data<-read.table("./data/log2_rslts_cph_hiv.pos_chf_geneName_corrected_noComplex.txt",header=TRUE,sep="\t")
data<-read.delim("../data/log2_rslts_cph_hiv.pos_mortality_geneName_corrected_noComplex.txt",header=TRUE,sep="\t")
dim(data)
colnames(data)



###some uniprot ids are not consistent. use symbol

colnames(data)[3]="pvalue"
colnames(data)[7]="Gene.Name"
NAs=which(is.na(data$Gene.Name)==TRUE)
if(length(NAs)>0){
  data=data[-NAs,]
}
data_sub=data%>%group_by(Gene.Name)%>%mutate('min.pvalue' = min(pvalue))%>%dplyr::select(min.pvalue,Gene.Name)
data_sub=unique(data_sub)

sum(table(data_sub$Gene.Name)>1)
head(sort((table(data_sub$Gene.Name)),decreasing = TRUE))
data_sub[which(data_sub$Gene.Name%in%names(head(sort((table(data_sub$Gene.Name)),decreasing = TRUE)))),]


data=as.data.frame(data)

colnames(data)[c(3,7)]=c("p.value","SYMBOL")
matched_entrez=match_entrez(data)
get_rsea(matched_entrez,go_list,go_annotation, output_name="../results/log2_rslts_cph_hiv.pos_mortality_geneName_corrected_noComplex_rsea_go")
get_rsea(matched_entrez,wp_list,wp_annotation, output_name="../results/log2_rslts_cph_hiv.pos_mortality_geneName_corrected_noComplex_rsea_wp")

jensen_annotation=jensen_annotation[,c(1,2,4)]
colnames(jensen_annotation)[3]="gene"
get_rsea(matched_entrez,jensen_list,jensen_annotation, output_name="../results/log2_rslts_cph_hiv.pos_mortality_geneName_corrected_noComplex_rsea_jensen")


#set sizes
hist(lengths(go_list), breaks=100)
hist(lengths(wp_list), breaks=100)
hist(lengths(jensen_list), breaks=100)
min(lengths(go_list))
max(lengths(go_list))
min(lengths(wp_list))
max(lengths(wp_list))
min(lengths(jensen_list))
max(lengths(jensen_list))

renal1_jensen<-readxl::read_excel("~/Dropbox (Gladstone)/GB-PH-1244/enrichment_results_v2/log2_rslts_cph_hiv.pos_renal1_geneName_corrected_noComplex_rsea_jensen_all_results.xlsx")
renal1_wp<-readxl::read_excel("~/Dropbox (Gladstone)/GB-PH-1244/enrichment_results_v2/log2_rslts_cph_hiv.pos_renal1_geneName_corrected_noComplex_rsea_wp_all_results.xlsx")
renal1_go<-readxl::read_excel("~/Dropbox (Gladstone)/GB-PH-1244/enrichment_results_v2/log2_rslts_cph_hiv.pos_renal1_geneName_corrected_noComplex_rsea_go_all_results.xlsx")

hist(renal1_go$Size, breaks=100)
hist(renal1_wp$Size, breaks=100)
hist(renal1_jensen$Size, breaks=100)
