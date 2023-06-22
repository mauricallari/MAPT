###		Main analyses in the TCGA	###

# WORKSPACE CLEARING
rm(list=ls()) 
graphics.off()

# LOAD REQUIRED LIBRARIES
library(affy)
library(reshape2)

# SET DIRECTORY
wd<-"/path/to/wd/"
setwd(wd)

# LOAD DEPMAP FILTERED DATASET
load(paste0(wd,"DEPMAP_gene_expression.RData"))
dataset<-eset


##########################################################
######	   	MAPT correlation with drug response    	######
##########################################################

# IMPORT DRUG RESPONSE DATA
data<-read.csv(paste0(wd,"AUC_per_cell_line_per_drug.csv"),header=T)
subdata<-data[,c("ARXSPAN_ID","DRUG_NAME","auc")]
mymean<-function(x){
  m<-mean(x,na.rm=T)
}
tab<-dcast(subdata,ARXSPAN_ID~DRUG_NAME,value.var="auc",fun.aggregate=mymean)
tab<-tab[-1,-2]
response.all<-tab[,2:ncol(tab)]
rownames(response.all)<-tab$ARXSPAN_ID

# cell lines common to GEP and DRUG data
common<-intersect(colnames(dataset),rownames(response.all))
length(common)
dataset.com<-dataset[,common]
response.com<-response.all[common,]

# CORRELATION MAPT VS RESPONSE
summary.linage<-summary(as.factor(dataset.com$lineage))
types<-names(summary.linage[summary.linage>=10])
cor.all<-numeric()
pval.all<-numeric()
filter.all<-logical()
for(i in 1:length(types)) {
 
  # subset
  subset<-dataset.com[,dataset.com$lineage==types[i]]
  sub.response<-response.com[colnames(subset),]
 
  # build drug filter
  sum.na<-apply(sub.response,2,function(x){sum(is.na(x)==F)})
  sum.sens<-apply(sub.response,2,function(x){sum(x<0.8,na.rm=T)})
  sum.res<-apply(sub.response,2,function(x){sum(x>0.8,na.rm=T)})
  filter<-ifelse(sum.na>=10&sum.sens>=2&sum.res>=2,1,NA)
  filter.all<-cbind(filter.all,filter)
  
  # correlation
  MAPT<-exprs(subset)[fData(subset)$Gene=="MAPT",]
  cor.type<-numeric()
  pval.type<-numeric()
  for(j in 1:ncol(sub.response)) {
    drug.resp<-sub.response[,j]
    if(sum(is.na(drug.resp))<(length(drug.resp)-2)) {
      cor<-cor.test(drug.resp,MAPT,method="spearman")$estimate
      pval<-cor.test(drug.resp,MAPT,method="spearman")$p.value
    } else {
      cor<-NA
      pval<-NA
    }
    cor.type<-c(cor.type,cor)
    pval.type<-c(pval.type,pval)
  }
  cor.all<-cbind(cor.all,cor.type)
  pval.all<-cbind(pval.all,pval.type)
}

colnames(cor.all)<-types
rownames(cor.all)<-colnames(response.com)
colnames(pval.all)<-types
rownames(pval.all)<-colnames(response.com)

cor.f<-cor.all*filter.all
pval.f<-pval.all*filter.all
sum.good<-apply(cor.f,1,function(x){sum(is.na(x)==F)})
sum.sig<-apply(pval.f,1,function(x){sum(x<0.05,na.rm=T)})
cor.good<-cor.f[sum.good>11&sum.sig>0,]
cor.good<-cor.good[rownames(cor.good)!="BLEOMYCIN (50 UM)",]
pval.good<-pval.f[sum.good>11&sum.sig>0,]
pval.good<-pval.good[rownames(pval.good)!="BLEOMYCIN (50 UM)",]

cor.output<-data.frame(Drug=rownames(cor.good),cor.good)
write.table(cor.output,file="Correlation_MAPT_response_pan_cancer.txt",quote=F,row.names=F,sep="\t")

pval.output<-data.frame(Drug=rownames(pval.good),pval.good)
write.table(pval.output,file="Pvalue_MAPT_response_pan_cancer.txt",quote=F,row.names=F,sep="\t")
