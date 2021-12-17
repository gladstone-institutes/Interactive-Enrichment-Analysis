
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



run_ORA5<-function(data, logFC_data, list_result, database, database_annotation){
  #data=as.matrix(pvalue_results_human_voom[,1]);list_result=ora_results_human_voom_wp;database=wp_list;database_annotation=wp_annotation
  
  
  GeneRatio=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  BgRatio=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  pvalue=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  p.adjust=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  qvalue=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  Count=matrix(NA,nrow=length(unique(database_annotation$set_id)),ncol=1)
  rownames(GeneRatio)=unique(database_annotation$set_id)
  rownames(BgRatio)=unique(database_annotation$set_id)
  rownames(pvalue)=unique(database_annotation$set_id)
  rownames(p.adjust)=unique(database_annotation$set_id)
  rownames(qvalue)=unique(database_annotation$set_id)
  rownames(Count)=unique(database_annotation$set_id)
  
  na_row=which(is.na(data)==TRUE)
  if(length(na_row)>0){
    data=as.matrix(data[-na_row,1]) 
  }
  
  data_m=as.data.frame(cbind(rownames(data),data))
  colnames(data_m)=c("Gene","pvalue")
  merged=merge(data_m,gene_entrez,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
  
  logFC_data=as.matrix(logFC_data)
  logFC_data=cbind(logFC_data,rownames(logFC_data))
  colnames(logFC_data)=c("logFC","Gene")
  
  merged=merge(merged,logFC_data,by="Gene",all = "L",all.x=TRUE,all.y=FALSE)
  merged=merged[!(is.na(merged$ENTREZID)),]
  enrichment_result <- clusterProfiler::enricher(
    merged$ENTREZID[ intersect( which(as.numeric(as.character(merged$pvalue)) < 0.05),  which( abs(as.numeric(as.character(merged$logFC))) > 1)) ],
    universe = merged$ENTREZID,
    pAdjustMethod = "holm",
    pvalueCutoff = 1, #p.adjust cutoff
    qvalueCutoff = 1,
    minGSSize = 1,
    maxGSSize = 100000,
    TERM2GENE = database_annotation[,c("set_id","gene")],
    TERM2NAME = database_annotation[,c("set_id","set_id")])
  
  #enrichment_result <- DOSE::setReadable(enrichment_result, org.Hs.eg.db, keyType = "ENTREZID")
  if(length(enrichment_result)>0){
    enrichment_result=as.data.frame(enrichment_result)
    enrichment_result=enrichment_result[,-match("geneID",colnames(enrichment_result))]
    
    GeneRatio[match(enrichment_result$ID,rownames(GeneRatio))]=enrichment_result$GeneRatio
    BgRatio[match(enrichment_result$ID,rownames(GeneRatio))]=enrichment_result$BgRatio
    pvalue[match(enrichment_result$ID,rownames(GeneRatio))]=enrichment_result$pvalue
    p.adjust[match(enrichment_result$ID,rownames(GeneRatio))]=enrichment_result$p.adjust
    qvalue[match(enrichment_result$ID,rownames(GeneRatio))]=enrichment_result$qvalue
    Count[match(enrichment_result$ID,rownames(GeneRatio))]=enrichment_result$Count
    
    list_result=modifyList(list_result, list(GeneRatio = cbind(list_result$GeneRatio,GeneRatio), BgRatio = cbind(list_result$BgRatio,BgRatio) ,
                                             pvalue = cbind(list_result$pvalue,pvalue), p.adjust = cbind(list_result$p.adjust,p.adjust),  
                                             qvalue = cbind(list_result$qvalue,qvalue), Count = cbind(list_result$Count,Count)) )
    
  }else{
    list_result=modifyList(list_result, list(GeneRatio = cbind(list_result$GeneRatio,GeneRatio), BgRatio = cbind(list_result$BgRatio,BgRatio) ,
                                             pvalue = cbind(list_result$pvalue,pvalue), p.adjust = cbind(list_result$p.adjust,p.adjust),  
                                             qvalue = cbind(list_result$qvalue,qvalue), Count = cbind(list_result$Count,Count)) )
    
  }
  
  return(list_result)
  
  
  
  
}







run_ORA_noLogFold<-function(data,  database, database_annotation){

  data=data[!(is.na(data$SYMBOL)),]
  data=data[!(is.na(data$ENTREZID)),]
  
  enrichment_result <- clusterProfiler::enricher(
    data$ENTREZID[  which(as.numeric(as.character(data$p.value)) < 0.05) ],
    universe = data$ENTREZID,
    #pAdjustMethod = "holm",
    pvalueCutoff = 1, #p.adjust cutoff
    qvalueCutoff = 1,
    minGSSize = 3,
    maxGSSize = 100000,
    TERM2GENE = database_annotation[,c("set_id","gene")],
    TERM2NAME = database_annotation[,c("set_id","set_id")])
  
 
  enrichment_result=as.data.frame(enrichment_result)
  
  result=NULL
 
    for(j in 1:nrow(enrichment_result)){
      annotation=database_annotation[which(database_annotation$set_id==enrichment_result[j,2])[1],1]
      all_genes=database_annotation[which(database_annotation$set_id==enrichment_result[j,2]),3]
      all_genes=data$SYMBOL[which(as.character(data$ENTREZID)%in%all_genes)]
   
      all_genes=paste(all_genes,collapse="/")
      
      gene_hits=database_annotation[which(database_annotation$set_id==enrichment_result[j,2]),3]
      gene_hits=data$SYMBOL[intersect(which(as.character(data$ENTREZID)%in%gene_hits),which(data$p.value<0.05))]
      
      gene_hits=paste(gene_hits,collapse="/")
      
      result=rbind(result,cbind(enrichment_result[j,],annotation,all_genes,gene_hits))
      
    }
    
  
    result=result[!is.na(result$pvalue),]
    
    
  return(result)
  
  
  
  
}









setwd("/Users/mingyoungshin/Dropbox (Gladstone)/BC-ZY-1231")


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

##############renal
data_renal<-read.delim("../data/log2_rslts_cph_hiv.pos_renal1_geneName_corrected_noComplex.txt",header=TRUE,sep="\t")
dim(data_renal)
colnames(data_renal)



###some uniprot ids are not consistent. use symbol

colnames(data_renal)[3]="pvalue"
colnames(data_renal)[7]="Gene.Name"
NAs=which(is.na(data_renal$Gene.Name)==TRUE)
if(length(NAs)>0){
  data_renal=data_renal[-NAs,]
}
data_renal_sub=data_renal%>%group_by(Gene.Name)%>%mutate('min.pvalue' = min(pvalue))%>%dplyr::select(min.pvalue,Gene.Name)
data_renal_sub=unique(data_renal_sub)

sum(table(data_renal_sub$Gene.Name)>1)
head(sort((table(data_renal_sub$Gene.Name)),decreasing = TRUE))
data_renal_sub[which(data_renal_sub$Gene.Name%in%names(head(sort((table(data_renal_sub$Gene.Name)),decreasing = TRUE)))),]


##############diabetes
data_diabetes<-read.delim("../data/log2_rslts_cph_hiv.pos_diabetes_geneName_corrected_noComplex.txt",header=TRUE,sep="\t")
dim(data_diabetes)
colnames(data_diabetes)



###some uniprot ids are not consistent. use symbol

colnames(data_diabetes)[3]="pvalue"
colnames(data_diabetes)[7]="Gene.Name"
NAs=which(is.na(data_diabetes$Gene.Name)==TRUE)
if(length(NAs)>0){
  data_diabetes=data_diabetes[-NAs,]
}
data_diabetes_sub=data_diabetes%>%group_by(Gene.Name)%>%mutate('min.pvalue' = min(pvalue))%>%dplyr::select(min.pvalue,Gene.Name)
data_diabetes_sub=unique(data_diabetes_sub)

sum(table(data_diabetes_sub$Gene.Name)>1)
head(sort((table(data_diabetes_sub$Gene.Name)),decreasing = TRUE))
data_diabetes_sub[which(data_diabetes_sub$Gene.Name%in%names(head(sort((table(data_diabetes_sub$Gene.Name)),decreasing = TRUE)))),]




##############chf
data_chf<-read.delim("../data/log2_rslts_cph_hiv.pos_chf_geneName_corrected_noComplex.txt",header=TRUE,sep="\t")
dim(data_chf)
colnames(data_chf)



###some uniprot ids are not consistent. use symbol

colnames(data_chf)[3]="pvalue"
colnames(data_chf)[7]="Gene.Name"
NAs=which(is.na(data_chf$Gene.Name)==TRUE)
if(length(NAs)>0){
  data_chf=data_chf[-NAs,]
}
data_chf_sub=data_chf%>%group_by(Gene.Name)%>%mutate('min.pvalue' = min(pvalue))%>%dplyr::select(min.pvalue,Gene.Name)
data_chf_sub=unique(data_chf_sub)

sum(table(data_chf_sub$Gene.Name)>1)
head(sort((table(data_chf_sub$Gene.Name)),decreasing = TRUE))
data_chf_sub[which(data_chf_sub$Gene.Name%in%names(head(sort((table(data_chf_sub$Gene.Name)),decreasing = TRUE)))),]


all_genes=unique(c(data_renal$UniProt,data_diabetes$UniProt,data_chf$UniProt))
sig_genes=unique(intersect(data_renal$UniProt[data_renal$q.value<0.05],data_diabetes$UniProt[data_diabetes$q.value<0.05]))
sig_genes=unique(intersect(sig_genes,data_chf$UniProt[data_chf$q.value<0.05]))
write.table(sig_genes,file="genes_significant_across3sets.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)


sig_genes2=unique(intersect(data_renal$UniProt[data_renal$pvalue<0.05],data_diabetes$UniProt[data_diabetes$pvalue<0.05]))
sig_genes2=unique(intersect(sig_genes2,data_chf$UniProt[data_chf$pvalue<0.05]))


bald_down_go_enrichment=run_ORA_noLogFold(bald_down,go_list,go_annotation)
jensen_list,jensen_annotation
