rm(list = ls())
myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");myfilepath
setwd(myfilepath)

rm(list = ls())
load("TCGA-BLCA.Rdata")
data<-cbind(Gene=rownames(TCGA_exp),TCGA_exp)
write.table(data,"TCGA-BLCA.xls",quote = F,sep = "\t",row.names = F)

rm(list = ls())
BLCA_Immune_activity<-read.table("BLCA- Immune.activity .txt",sep = "\t",header = T,check.names = F) 
BLCA_Immune_activity1<-BLCA_Immune_activity[-c(1:3,21:23),]
BLCA_Immune_activity1$Steps<-gsub("Step4.","",BLCA_Immune_activity1$Steps)
rownames(BLCA_Immune_activity1)<-BLCA_Immune_activity1$Steps
BLCA_Immune_activity2<-BLCA_Immune_activity1[,-1]
for (i in 1:ncol(BLCA_Immune_activity2)) {
  #i=1
  a<-unlist(strsplit(names(BLCA_Immune_activity2)[i],""));a
  names(BLCA_Immune_activity2)[i]<-paste(a[1:16],collapse = "")
}
names(BLCA_Immune_activity2)
BLCA_Immune_activity3<-data.frame(t(BLCA_Immune_activity2))
rownames(BLCA_Immune_activity3)<-gsub("\\.","-",rownames(BLCA_Immune_activity3))

BLCA_Immune_activity3$GSM<-rownames(BLCA_Immune_activity3)

load("TCGA_KM_data.Rdata")
nam<-intersect(BLCA_Immune_activity3$GSM,KM_data$GSM)
KM_data<-KM_data[,c(1,5)]
data<-merge(BLCA_Immune_activity3,KM_data,by="GSM")
library(reshape2)
TME_New<-melt(data[,-1])
head(TME_New)
colnames(TME_New)=c("Group","gene","Expression")  #设置行名
# 3.3 出图
library(ggthemes)
library(ggplot2)
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }
library(ggpubr)
head(TME_New)
# TME_New$Expression<-log2(TME_New$Expression+0.0001)
box_TME<-ggboxplot(TME_New, x="gene", y="Expression", color = "Group", 
                   palette=c("#f85a40","#00b2a9"),
                   #add = "jitter"
)+
  labs(y="Immune_activity",x= NULL,title = "")+  
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = F);box_TME
ggsave("TIP.pdf",box_TME,height=5,width=10)
