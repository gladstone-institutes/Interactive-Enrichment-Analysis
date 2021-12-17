
database_lists1<-load("/Users/mingyoungshin/Dropbox (Gladstone)/PFOCRInPathwayAnalyses_git/PFOCRInPathwayAnalyses/databases.RData")#has wp, pfocr, go
#database_lists1=c(database_lists1,load("/Users/mingyoungshin/Dropbox (Gladstone)/enrichment_databases/jensen_database.RData"))

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

get_gsea<-function(score_calculated = FALSE, data, database, database_annotation, output_name="temp"){
  #score_calculated = TRUE; data=matched_entrez;dabase= go_list;database_annotation=go_annotation;output_name="../HFrEF_gsea_go"
  ### geneList prep
  if(sorted_gene_list==FALSE){
    
    gene_list <- data %>%
      mutate(Score = (-1*(sign(as.numeric(log2FC)) ))* - log10(as.numeric(as.character(p.value)))) %>%
      dplyr::select(c("Score","ENTREZID")) %>%
      drop_na(ENTREZID) %>%
      dplyr::arrange(desc(Score))
    gene_list <- unlist(split(gene_list[, 1], gene_list[, 2]))
    
    gene_list = sort(gene_list[unique(names(gene_list))], decreasing = TRUE)
    
  }else{
    
    gene_list <- data %>%
      dplyr::select(c("Score","ENTREZID")) %>%
      drop_na(ENTREZID) %>%  drop_na(Score) %>%
      dplyr::arrange(desc(Score))
    gene_list <- unlist(split(gene_list[, 1], gene_list[, 2]))
    
    gene_list = sort(gene_list[unique(names(gene_list))], decreasing = TRUE)
  
  }
  #data=up
  
  
  
  enrichment_result <- GSEA(
    gene_list,
  #  pAdjustMethod="holm",
    TERM2GENE = database_annotation[,c("set_id","gene")],
    TERM2NAME = database_annotation[,c("set_id","set_id")]    ,
    minGSSize = 3,
    maxGSSize = 100000,
    pvalueCutoff = 1,
    verbose=FALSE)
 
  
  enrichment_result=as.data.frame(enrichment_result)
  
  core_enrichment=sapply(enrichment_result$core_enrichment,function(x){
    x=strsplit(as.character(x), '/')
    paste(data[match(unlist(x),as.character(data[,"ENTREZID"])),"gene"],collapse = "/")
  })
  length(core_enrichment)
  dim(enrichment_result)
  enrichment_result$core_enrichment=core_enrichment
  
  result=NULL
  
  for(j in 1:nrow(enrichment_result)){
    annotation=database_annotation[which(database_annotation$set_id==enrichment_result[j,2])[1],1]
    genes=database_annotation[which(database_annotation$set_id==enrichment_result[j,2]),3]
    genes=data$gene[which(as.character(data$ENTREZID)%in%as.character(genes))]
    
    all_genes=paste(genes,collapse="/")
    
    
    
    if(score_calculated==FALSE){
      gene_hits=database_annotation[which(database_annotation$set_id==enrichment_result[j,2]),3]
      gene_hits=data$gene[intersect(which(as.character(data$ENTREZID)%in%as.character(gene_hits)),which(data$pvalue<0.05))]
      
      gene_hits=paste(gene_hits,collapse="/")
      result=rbind(result,cbind(enrichment_result[j,],annotation,all_genes,gene_hits))
    }else{
      result=rbind(result,cbind(enrichment_result[j,],annotation,all_genes))
    }
 
    
    
    
  }
  
  result=result[!is.na(result$pvalue),]
  result=result[order(result$pvalue,decreasing = FALSE),]
 
  write.table(result,paste0(output_name,"_all_results.txt"),col.names=colnames(result),row.names=FALSE,sep="\t",quote=FALSE)
  
  return(result)
}

match_entrez<-function(data, id_type="SYMBOL"){
  gene_entrez <- clusterProfiler::bitr(data$gene, fromType = id_type,
                                       toType = c("ENTREZID"),
                                       OrgDb = org.Hs.eg.db, drop=FALSE)
  
  tb=table(gene_entrez[,id_type])
  print(max(tb[which(tb>1)]))
  
  if(max(tb[which(tb>1)])>1){
    
    remove_i=NULL
    for(n in names(tb[which(tb>1)])){
      all_i=which(n==gene_entrez[,id_type])
      remove_i=c(remove_i,all_i[-1])
    }
    
  }
  print(dim(gene_entrez))
  if(length(remove_i)>0) gene_entrez=gene_entrez[-remove_i,]
  
  
  colnames(data)[1]=id_type
  print(dim(gene_entrez))
  print(dim(data))
  data=merge(data,gene_entrez,by=id_type)
  print(dim(data))
  colnames(data)[ncol(data)]="ENTREZID"
  
  print(dim(data))
  return(data)
}

match_symbol<-function(data, id_type="UNIPROT"){
  gene_symbol <- clusterProfiler::bitr(data[,id_type], fromType = id_type,
                                       toType = c("SYMBOL"),
                                       OrgDb = org.Hs.eg.db, drop=FALSE)
  
  tb=table(gene_symbol[,id_type])
  print(max(tb[which(tb>1)]))
  
  if(max(tb[which(tb>1)])>1){

    remove_i=NULL
    for(n in names(tb[which(tb>1)])){
      all_i=which(n==gene_symbol[,id_type])
      remove_i=c(remove_i,all_i[-1])
    }

  }
  print(dim(gene_symbol))
  if(length(remove_i)>0) gene_symbol=gene_symbol[-remove_i,]
  

  colnames(data)[1]=id_type
  print(dim(gene_symbol))
  print(dim(data))
  data=merge(data,gene_symbol,by=id_type)
  print(dim(data))
  colnames(data)[ncol(data)]="SYMBOL"
  
  print(dim(data))
  return(data)
}





library(dplyr)
library(tidyr)
library(magrittr)
library("org.Hs.eg.db")
library(writexl)
library(clusterProfiler)


data<-read.csv("/Users/mingyoungshin/Dropbox (Gladstone)/20210701_Dubin_CDK_Proteomics/Enrichment/HFrEF_ranked.csv",header=TRUE)
dim(data)
colnames(data)

colnames(data)=c("gene","Score")

matched_entrez=match_entrez(data, "UNIPROT")
matched_entrez=matched_entrez[!is.na(matched_entrez$ENTREZID),]# 0.76% of input gene IDs are fail to map
dim(matched_entrez)
head(matched_entrez)

colnames(matched_entrez)[1]="gene"
HFrEF_wp=get_gsea(score_calculated = TRUE, matched_entrez, wp_list,wp_annotation, output_name="./HFrEF_gsea_wp")
HFrEF_pfocr=get_gsea(score_calculated = TRUE, matched_entrez, pfocr_list,pfocr_annotation, output_name="./HFrEF_gsea_pfocr")
HFrEF_go=get_gsea(score_calculated = TRUE, matched_entrez, go_list,go_annotation, output_name="./HFrEF_gsea_go")

write_xlsx(HFrEF_go,"/HFrEF_gsea_go__all_results.xlsx")






data<-read.csv("/Users/mingyoungshin/Dropbox (Gladstone)/20210701_Dubin_CDK_Proteomics/Enrichment/generalHF_ranked.csv",header=TRUE)
dim(data)
colnames(data)

colnames(data)=c("gene","Score")

matched_entrez=match_entrez(data, "UNIPROT")
matched_entrez=matched_entrez[!is.na(matched_entrez$ENTREZID),]# 0.76% of input gene IDs are fail to map
dim(matched_entrez)
head(matched_entrez)

colnames(matched_entrez)[1]="gene"
generalHF_wp=get_gsea(score_calculated = TRUE, matched_entrez, wp_list,wp_annotation, output_name="./generalHF_gsea_wp")
generalHF_pfocr=get_gsea(score_calculated = TRUE, matched_entrez, pfocr_list,pfocr_annotation, output_name="./generalHF_gsea_pfocr")
generalHF_go=get_gsea(score_calculated = TRUE, matched_entrez, go_list,go_annotation, output_name="./generalHF_gsea_go")






data<-read.csv("/Users/mingyoungshin/Dropbox (Gladstone)/20210701_Dubin_CDK_Proteomics/Enrichment/HFpEF_ranked.csv",header=TRUE)
dim(data)
colnames(data)

colnames(data)=c("gene","Score")

matched_entrez=match_entrez(data, "UNIPROT")
matched_entrez=matched_entrez[!is.na(matched_entrez$ENTREZID),]# 0.76% of input gene IDs are fail to map
dim(matched_entrez)
head(matched_entrez)

colnames(matched_entrez)[1]="gene"
HFpEF_wp=get_gsea(score_calculated = TRUE, matched_entrez, wp_list,wp_annotation, output_name="./HFpEF_gsea_wp")
HFpEF_pfocr=get_gsea(score_calculated = TRUE, matched_entrez, pfocr_list,pfocr_annotation, output_name="./HFpEF_gsea_pfocr")
HFpEF_go=get_gsea(score_calculated = TRUE, matched_entrez, go_list,go_annotation, output_name="./HFpEF_gsea_go")

