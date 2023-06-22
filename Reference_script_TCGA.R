###		Main analyses in the TCGA	###

# WORKSPACE CLEARING
rm(list=ls()) 
graphics.off()

# LOAD REQUIRED LIBRARIES
library(affy)
library(phenoTest)

# SET DIRECTORY
wd<-"/path/to/wd/"
setwd(wd)

# LOAD TCGA FILTERED DATASET
load(paste0(wd,"TCGA_pancancer_RNAseq_mutAnnotated_filtered.RData"))
dataset<-gep.f


##########################################################
######	    	MAPT correlation with all genes      	######
##########################################################
types<-unique(dataset$acronym)
cor.all<-numeric()
for(i in 1:length(types)) {
  subset<-dataset[,which(dataset$acronym==types[i])]
  MAPT<-exprs(subset)[fData(subset)$Symbol=="MAPT",]
  others<-exprs(subset)[fData(subset)$Symbol!="MAPT",]
  cor<-apply(others,1,function(x){cor(x,MAPT,method="spearman")})  
  cor.all<-cbind(cor.all,cor)
  colnames(cor.all)[i]<-paste0(types[i])
}

output<-data.frame(Gene=rownames(cor.all),cor.all)
output[is.na(output)]<-0
write.table(output,file="MAPT_AllGenes_correlation_pancancer.txt",row.names=F,sep="\t",quote=F)
##########################################################


##########################################################
######	     	  GENESET ENRICHMENT ANALYSIS        	######
##########################################################
# select one gene symbol 
cor.all<-cor.all[grep("\\?",rownames(cor.all),invert=T),]
symbols<-apply(as.matrix(rownames(cor.all)),1,function(x){strsplit(x,"\\|")[[1]][1]})
rownames(cor.all)<-symbols

# load genesets
gsets.df<-read.table(paste0(wd,"H.all_v7.1_panther_custom.txt"), sep="\t",
                     fill = TRUE, stringsAsFactors = FALSE, header = F)

# shorten long names, filter and and create geneset list
mysub<-function(x) {
  if(nchar(x)>100) {x<-substring(x,1,100)}
  return(x)
}
gsets.names<-apply(as.matrix(gsets.df[, 1]),1,mysub)
gsets.names<-gsub("HALLMARK_","H_",gsets.names)
gsets.names<-gsub("PANTHER_","P_",gsets.names)
gsets.genes<-gsets.df[, -c(1)]
gsets.list<-split(gsets.genes, seq_len(nrow(gsets.genes)))
gsets.list<-lapply(gsets.list, function(x) x[x != ""])
names(gsets.list)<-gsets.names
gsets.length<-sapply(gsets.list,"length")
filter<-gsets.length>=10&gsets.length<=200
gsets.list<-gsets.list[filter]

# GSEA
set.seed(123)
nes.all<-numeric()
fdr.all<-numeric()
for(j in 1:ncol(cor.all)) {
  # run gsea and export summary
  x<-cor.all[,j]
  names(x)<-rownames(cor.all)
  gsea<-gsea(x,gsets=gsets.list,
             logScale=F,p.adjust.method="BH",pval.comp.method="signed",
             minGenes=0,maxGenes=1000,center=T)
  res<-summary.gseaData(gsea)
  nes.all<-cbind(nes.all,res[,"nes"])
  fdr.all<-cbind(fdr.all,res[,"fdr"])
}  
colnames(fdr.all)<-colnames(cor.all) 
colnames(nes.all)<-colnames(cor.all) 

fdr.output<-data.frame(Geneset=rownames(fdr.all),fdr.all)  
write.table(fdr.output,file="GSEA_FDR.txt",quote=F,sep="\t",row.names=F)

nes.output<-data.frame(Geneset=rownames(nes.all),nes.all)  
write.table(nes.output,file="GSEA_NES.txt",quote=F,sep="\t",row.names=F)
##########################################################


##########################################################
######	     	  SURVIVAL ANALYSIS        	######
##########################################################
# subset data with survival info
subset<-dataset[,dataset$vital_status%in%c("Alive","Dead")]

# combine P53 mutations and CNA
subset$mut_cna<-paste(subset$TP53_MUTATIONS,subset$TP53_CNA, sep=";")
subset$p53_mut<-F
subset$p53_mut[subset$mut_cna%in%c(";homdel_rec",
                                   "Inframe Mutation (putative driver);",
                                   "Inframe Mutation (putative driver);homdel_rec",
                                   "Missense Mutation (putative driver);",
                                   "Missense Mutation (putative driver);Amplification",
                                   "Missense Mutation (putative driver);homdel_rec",
                                   "Truncating mutation (putative driver);",
                                   "Truncating mutation (putative driver);Amplification",
                                   "Truncating mutation (putative driver);homdel_rec")]<-T
subset$p53_mut[is.na(subset$TP53_MUTATIONS)]<-NA

# create classes for clinico-path variables 
Size<-subset$pathologic_T
Size[Size%in%c("Tis","T0","T1","T1a","T1a1","T1b","T1b1","T1b2","T1c","T2","T2a","T2a1","T2a2","T2b","T2c")]<-"T0T2"
Size[Size%in%c("T3","T3a","T3b","T3c","T4","T4a","T4b","T4c","T4d","T4e")]<-"T3T4"
Size[Size%in%c("","[Discrepancy]","[Not Applicable]","[Not Available]","TX")]<-NA

LN<-subset$pathologic_N
LN[LN%in%c("N0","N0 (i-)","N0 (i+)","N0 (mol+)")]<-"N0"
LN[LN%in%c("N1","N1a","N1b","N1c","N1mi","N2","N2a","N2b","N2c","N3","N3a","N3b","N3c")]<-"N1"
LN[LN%in%c("","[Discrepancy]","[Not Applicable]","[Not Available]","NX")]<-NA

Met<-subset$pathologic_M
Met[Met%in%c("cM0 (i+)","M0")]<-"M0"
Met[Met%in%c("M1","M1a","M1b","M1c")]<-"M1"
Met[Met%in%c("","[Discrepancy]","[Not Applicable]","[Not Available]","MX","[Unknown]")]<-NA

cancer<-subset$acronym

event<-ifelse(subset$vital_status=="Dead",1,0)
OS<-subset$days_to_death
OS[event==0]<-subset$days_to_last_followup[event==0]
OS<-as.numeric(OS)/365.25
OS[which(OS>=5)]<-5
event[which(OS>=5)]<-0

MAPT<-exprs(subset)[fData(subset)$Symbol=="MAPT",]
AURKA<-exprs(subset)[fData(subset)$Symbol=="AURKA",]
p53_mut<-subset$p53_mut

# suvival analysis dataframe
df<-data.frame(cancer,p53_mut,event,OS,MAPT,Size,LN,Met,AURKA)
df<-df[is.na(OS)==F&OS>0,]

# COX ANALYSIS MAPT AND OS
cancer.list<-unique(as.vector(df$cancer))
survival.all<-matrix(nrow=length(cancer.list),ncol=18)
for(i in 1:length(cancer.list)) {
  # all samples univar
  subdf<-df[df$cancer==cancer.list[i],]
  med.MAPT<-median(subdf$MAPT,na.rm=T)
  subdf$MAPTcat<-ifelse(subdf$MAPT>med.MAPT,"H","L")
  subdf$MAPTcat<-factor(subdf$MAPTcat,levels=c("L","H"))
  if(sum(subdf$event)>=10) {
    cox.uni<-coxph(Surv(OS,event)~MAPTcat,data=subdf)
    survival.all[i,1]<-cox.uni$n
    survival.all[i,2]<-summary(cox.uni)$coefficients[1,2]
    survival.all[i,3]<-summary(cox.uni)$coefficients[1,5]
 
    # all samples multivar
    var.unique<-sapply(lapply(subdf[,6:9], unique), function(x){length(levels(as.factor(unique(x))))})
    var.notNA<-apply(subdf[,6:9], 2, function(x){sum(is.na(x)==F)})
    var.include<-names(var.unique)[var.unique>1&var.notNA>20]
    if(length(var.include>=3)) {
      model.string<-paste0("Surv(OS,event)~MAPTcat+",paste(var.include,collapse = "+"))
      cox.multi<-coxph(as.formula(model.string),data=subdf)
      survival.all[i,4]<-cox.multi$n
      survival.all[i,5]<-summary(cox.multi)$coefficients["MAPTcatH",2]
      survival.all[i,6]<-summary(cox.multi)$coefficients["MAPTcatH",5]
    } 
  }
  
  # stratified by TP53 status
  if(sum(subdf$p53_mut,na.rm=T)>20 & sum(subdf$p53_mut==F,na.rm=T)>20) {
    # WT
    subdf.wt<-df[which(df$cancer==cancer.list[i]&df$p53_mut==F),]
    med.wt.MAPT<-median(subdf.wt$MAPT,na.rm=T)
    subdf.wt$MAPTcat<-ifelse(subdf.wt$MAPT>med.wt.MAPT,"H","L")
    subdf.wt$MAPTcat<-factor(subdf.wt$MAPTcat,levels=c("L","H"))   
    if(sum(subdf.wt$event)>=10) {
      cox.uni.wt<-coxph(Surv(OS,event)~MAPTcat,data=subdf.wt)
      survival.all[i,7]<-cox.uni.wt$n
      survival.all[i,8]<-summary(cox.uni.wt)$coefficients[1,2]
      survival.all[i,9]<-summary(cox.uni.wt)$coefficients[1,5]
      
      # WT samples multivar
      var.unique<-sapply(lapply(subdf.wt[,6:9], unique), function(x){length(levels(as.factor(unique(x))))})
      var.notNA<-apply(subdf.wt[,6:9], 2, function(x){sum(is.na(x)==F)})
      var.include<-names(var.unique)[var.unique>1&var.notNA>20]
      if(length(var.include>=3)) {
        model.string<-paste0("Surv(OS,event)~MAPTcat+",paste(var.include,collapse = "+"))
        cox.multi<-coxph(as.formula(model.string),data=subdf.wt)
        survival.all[i,10]<-cox.multi$n
        survival.all[i,11]<-summary(cox.multi)$coefficients["MAPTcatH",2]
        survival.all[i,12]<-summary(cox.multi)$coefficients["MAPTcatH",5]
      } 
    }
    
    # MUT
    subdf.mut<-df[which(df$cancer==cancer.list[i]&df$p53_mut==T),]
    med.mut.MAPT<-median(subdf.mut$MAPT,na.rm=T)
    subdf.mut$MAPTcat<-ifelse(subdf.mut$MAPT>med.mut.MAPT,"H","L")
    subdf.mut$MAPTcat<-factor(subdf.mut$MAPTcat,levels=c("L","H"))      
    if(sum(subdf.mut$event)>=10) {
      cox.uni.mut<-coxph(Surv(OS,event)~MAPTcat,data=subdf.mut)
      survival.all[i,13]<-cox.uni.mut$n
      survival.all[i,14]<-summary(cox.uni.mut)$coefficients[1,2]
      survival.all[i,15]<-summary(cox.uni.mut)$coefficients[1,5]
      
      # MUT samples multivar
      var.unique<-sapply(lapply(subdf.mut[,6:9], unique), function(x){length(levels(as.factor(unique(x))))})
      var.notNA<-apply(subdf.mut[,6:9], 2, function(x){sum(is.na(x)==F)})
      var.include<-names(var.unique)[var.unique>1&var.notNA>20]
      if(length(var.include>=3)) {
        model.string<-paste0("Surv(OS,event)~MAPTcat+",paste(var.include,collapse = "+"))
        cox.multi<-coxph(as.formula(model.string),data=subdf.mut)
        survival.all[i,16]<-cox.multi$n
        survival.all[i,17]<-summary(cox.multi)$coefficients["MAPTcatH",2]
        survival.all[i,18]<-summary(cox.multi)$coefficients["MAPTcatH",5]
      } 
    } 
  }
}

rownames(survival.all)<-cancer.list
colnames(survival.all)<-c("n.uni","uni.hr","uni.pval",
                          "n.multi","multi.hr","multi.pval",
                          "n.uni.wt","uni.hr.wt","uni.pval.wt",
                          "n.multi.wt","multi.hr.wt","multi.pval.wt",
                          "n.uni.mut","uni.hr.mut","uni.pval.mut",
                          "n.multi.mut","multi.hr.mut","multi.pval.mut")

surv.output<-data.frame(Cancer=rownames(survival.all),survival.all)
write.table(surv.output,file="MAPT_Survival_analysis_uni_multi_TP53status.txt",row.names=F,sep="\t",quote=F)
