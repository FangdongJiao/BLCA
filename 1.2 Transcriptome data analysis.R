rm(list = ls())
myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");myfilepath
setwd(myfilepath)

#TCGA-BLCA
{
  rm(list = ls())
  library(data.table)
  library(dplyr)
  stad.fkpm=fread("TCGA-BLCA.htseq_fpkm.tsv.gz",
                  header = T, sep = '\t',data.table = F)
  stad.pro=fread("gencode.v22.annotation.gene.probeMap",
                 header = T, sep = '\t',data.table = F)
  stad.pro=stad.pro[,c(1,2)]
  stad.fkpm.pro=merge(stad.pro,stad.fkpm,
                      by.y ="Ensembl_ID",
                      by.x = "id" )
  # rownames(stad.fkpm.pro)=stad.fkpm.pro$gene
  library(tidyverse)
  stad.fkpm.pro=distinct(stad.fkpm.pro,gene,.keep_all = T)
  stad.fkpm.pro <- column_to_rownames(stad.fkpm.pro,"gene")
  stad.fkpm.pro<-stad.fkpm.pro[,-1]
  boxplot(stad.fkpm.pro[1:50,1:50])
  
  TCGA_exp<-stad.fkpm.pro
  boxplot(TCGA_exp[1:50,1:50])
  
  colnames(stad.fkpm.pro)
  x<-substr(colnames(stad.fkpm.pro),14,15)
  x<-as.numeric(x)
  
  y<-ifelse(x < 10,'tumor','normal')
  y<-factor(y)
  group_list=y
  table(group_list)
  TCGA_group<-data.frame(group_list)
  TCGA_group$GSM<-colnames(stad.fkpm.pro)
  table(TCGA_group$group_list)
  tumor<-TCGA_group[TCGA_group$group_list=="tumor",]
  
  
  stad.phe=fread("TCGA-BLCA.GDC_phenotype.tsv.gz",
                 header = T, sep = '\t',data.table = F)
  rownames(stad.phe)<-stad.phe$submitter_id.samples
  pdata<-stad.phe[colnames(TCGA_exp),]
  
  #重新确定group分组
  group1<-data.frame(cbind(pdata$submitter_id.samples,pdata$sample_type.samples))
  table(group1[,2])
  
  phenotype<-cbind(pdata$submitter_id.samples,
                   pdata$age_at_initial_pathologic_diagnosis,
                   pdata$gender.demographic,
                   pdata$tobacco_smoking_history,
                   pdata$lymphovascular_invasion_present,
                   pdata$tumor_stage.diagnoses,
                   pdata$clinical_T,
                   pdata$neoplasm_histologic_grade,
                   pdata$primary_therapy_outcome_success
  )
  phenotype<-data.frame(phenotype)
  names(phenotype)<-c("GSM",
                      "age","gender","smoking",
                      "lymphatic_invasion",
                      "stage","clinical_T",
                      "grade","therapy_outcome" )
  data<-phenotype
  data[is.na(data)]<-""
  df1<-data
  table(df1$age)
  table(df1$stage)
  df1$stage[df1$stage=="not reported"]<-NA
  df1$stage<-gsub("stage ","",df1$stage)
  # df1$stage[df1$stage=="ia"]<-"i"
  # df1$stage[df1$stage=="iia"]<-"ii"
  # df1$stage[df1$stage=="iib"]<-"ii"
  # df1$stage[df1$stage=="iic"]<-"ii"
  # df1$stage[df1$stage=="iiia"]<-"iii"
  # df1$stage[df1$stage=="iiib"]<-"iii"
  # df1$stage[df1$stage=="iiic"]<-"iii"
  # df1$stage[df1$stage=="iva"]<-"iv"
  # df1$stage[df1$stage=="ivb"]<-"iv"
  table(df1$stage)
  
  
  metdata<-read.table("TCGA-BLCA.survival.tsv",header = T)
  rownames(metdata)<-metdata$sample
  y<-tumor$GSM
  survival<-metdata[y,]
  survival<-survival[,-3]
  survival<-na.omit(survival)
  
  load("protein_coding.gtf.Rda")
  TCGA_exp<-TCGA_exp[protein_coding.gtf$gene_name,]
  TCGA_exp<-na.omit(TCGA_exp)
  
  save(TCGA_exp,TCGA_group,survival,df1,file='TCGA-BLCA.Rdata')
  load("TCGA-BLCA.Rdata")

}
