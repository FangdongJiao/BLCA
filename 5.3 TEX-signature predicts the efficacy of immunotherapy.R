rm(list = ls())
myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");myfilepath
setwd(myfilepath)

rm(list = ls())
load("TCGA_KM_data.Rdata")

#散点图
load("TCGA-BLCA.Rdata")
table(df1$therapy_outcome)
df1$therapy_outcome[df1$therapy_outcome=="Complete Remission/Response"]<-"CR/PR"
df1$therapy_outcome[df1$therapy_outcome=="Partial Remission/Response"]<-"CR/PR"
df1$therapy_outcome[df1$therapy_outcome=="Progressive Disease"]<-"PD/SD"
df1$therapy_outcome[df1$therapy_outcome=="Stable Disease"]<-"PD/SD"
df1<-df1[,c(1,9)]
df1[df1==""]<-NA
df1<-na.omit(df1)

fg<-merge(KM_data,df1,by="GSM")
fgg<-fg[,c(1,4,5,6)]
library(reshape2)
num<-melt(fgg)
names(num)
head(num)
library(ggbeeswarm)
library(ggpubr)
p<-ggplot(num,aes(x=therapy_outcome,y=value))+geom_boxplot()+
  geom_beeswarm(aes(color=factor(therapy_outcome)),cex=2,alpha = 0.5)+
  theme_bw()+ xlab("")+ylab("Risk score")+
  #scale_color_manual(values=c("red","blue"))+
  guides(color=F)+
  stat_compare_means(aes(group =  therapy_outcome),
                     label = "p.signif",
                     hide.ns = T,
                     label.x = 1.5);p
dev.off()
pdf("CPRPPDSD_point.pdf",width=2,height=3,onefile = F);p;dev.off()

# roc
library(pROC)
head(fgg)
roc1<-roc(fgg$therapy_outcome, fgg$score)
plot(smooth(roc1),col="red",
     legacy.axes=T,print.auc=TRUE)
dev.off()
pdf("CPRPPDSD_AUC.pdf",width =4,height = 4,onefile = F)
plot(smooth(roc1),col="blue")
legend("bottomright",legend="Risk score AUC:0.70",lty=1,  bty = "n")
dev.off()
