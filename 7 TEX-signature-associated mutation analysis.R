rm(list = ls())
myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");myfilepath
setwd(myfilepath)

rm(list = ls())
load("TCGA-BLCA_SNP.Rdata")
df<-data
rm(data)
sample<-data.frame(df$Tumor_Sample_Barcode,
                   df$Hugo_Symbol)
for (i in 1:length(sample$df.Tumor_Sample_Barcode)) {
  a<-sample$df.Tumor_Sample_Barcode[i]
  a1<-unlist(strsplit(a,""))
  a2<-a1[1:16]
  a3<-paste(a2, sep = "",collapse="")
  sample$GSM[i]<-a3
}  
df$GSM<-sample$GSM
laml.maf1<-df

load("TCGA_KM_data.Rdata")
KM_data<-KM_data[,c(1,4,5)]
high<-KM_data[KM_data$levle=="high",]
low<-KM_data[KM_data$levle=="low",]
library(dplyr)
laml.maf_low<-merge(laml.maf1,low,by="GSM")
laml.maf_high<-merge(laml.maf1,high,by="GSM")

#https://mp.weixin.qq.com/s/WG4JHs9RSm5IEJiiGEzDkg
#http://www.360doc.com/content/20/1201/22/65172408_948993557.shtml

#L临床数据
library(TCGAbiolinks)
library(tidyverse)
paad.clin <- GDCquery_clinic("TCGA-BLCA")
clin <- paad.clin %>% 
  dplyr::select(c(
    "submitter_id", "days_to_last_follow_up", "vital_status",
    "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m",
    "ajcc_pathologic_stage", "gender"
  )) %>%
  `names<-`(c("Tumor_Sample_Barcode", "time", "status", "T", "N", "M", "stage", "gender"))
clin$status[clin$status=="Alive"]<-0
clin$status[clin$status=="Dead"]<-1
clin$status[clin$status=="Not Reported"]<-NA
save(laml.maf_low,laml.maf_high,clin,file="tubianCNV.Rdata")

rm(list = ls())
load("tubianCNV.Rdata")
#high low 生存数据
for (i in 1:length(laml.maf_low$Tumor_Sample_Barcode)) {
  a<-laml.maf_low$Tumor_Sample_Barcode[i]
  a1<-unlist(strsplit(a,""))
  a2<-a1[1:12]
  a3<-paste(a2, sep = "",collapse="")
  laml.maf_low$Tumor_Sample_Barcode[i]<-a3
}

for (i in 1:length(laml.maf_high$Tumor_Sample_Barcode)) {
  a<-laml.maf_high$Tumor_Sample_Barcode[i]
  a1<-unlist(strsplit(a,""))
  a2<-a1[1:12]
  a3<-paste(a2, sep = "",collapse="")
  laml.maf_high$Tumor_Sample_Barcode[i]<-a3
}
names(clin)[1]
clin_low<-merge(laml.maf_low,clin,by="Tumor_Sample_Barcode")
laml.maf_low1<-clin_low
clin_low1<-data.frame(clin_low$Tumor_Sample_Barcode,
                      clin_low$stage,
                     clin_low$gender,
                     clin_low$status)
names(clin_low1)<-c("Tumor_Sample_Barcode",
                    "stage",
                    "gender","status")

clin_high<-merge(laml.maf_high,clin,by="Tumor_Sample_Barcode")
laml.maf_high1<-clin_high
clin_high1<-data.frame(clin_high$Tumor_Sample_Barcode,
                      clin_high$stage,
                      clin_high$gender,
                      clin_high$status)
names(clin_high1)<-c("Tumor_Sample_Barcode",
                     "stage",
                     "gender","status")
table(clin_high1$stage)
table(clin_high1$gender)
table(clin_high1$status)

#画MAF文件的summary图
library(maftools)
maf_high = read.maf(maf = laml.maf_high1,clinicalData = clin_high1)
getClinicalData(maf_high)

maf_low = read.maf(maf = laml.maf_low1,clinicalData = clin_low1)
getClinicalData(maf_low)
save(maf_high,maf_low,file="maf.Rdata")

rm(list = ls())
load("maf.Rdata")
#1\2\oncoplot for top 20 genes
vc_cols = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
pdf("low_oncoplot_20_genes.pdf",onefile = F,width = 10,height = 15)
oncoplot(maf = maf_low, titleText ="Low Risk Group",
         clinicalFeatures=c("stage","gender","status"),
         colors = vc_cols, top = 20,
         draw_titv = F)#
dev.off()
pdf("high_oncoplot_20_genes.pdf",onefile = F,width = 10,height = 15)
oncoplot(maf = maf_high, titleText ="High Risk Group",
         clinicalFeatures=c("stage","gender","status"),
         colors = vc_cols, top = 20,
         draw_titv = F)
dev.off()
