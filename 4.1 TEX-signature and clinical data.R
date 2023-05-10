myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");
setwd(myfilepath)

rm(list = ls())
load("TCGA_KM_data.Rdata")
table(KM_data$levle)
KM_data<-na.omit(KM_data)
KM_data<-KM_data[,c(1,4)]

load("TCGA-BLCA.Rdata")
data<-merge(df1,KM_data,by="GSM")
df<-data[,-c(1,4,7)]

df$age<-ifelse(df$age>60,">60","<=60")
df$therapy_outcome[df$therapy_outcome=="Complete Remission/Response"]<-"CR/PR"
df$therapy_outcome[df$therapy_outcome=="Partial Remission/Response"]<-"CR/PR"
df$therapy_outcome[df$therapy_outcome=="Progressive Disease"]<-"PD/SD"
df$therapy_outcome[df$therapy_outcome=="Stable Disease"]<-"PD/SD"

fff<-df
fff$stage[fff$stage=="i"]<-"I+II"
fff$stage[fff$stage=="ii"]<-"I+II"
fff$stage[fff$stage=="iii"]<-"III+IV"
fff$stage[fff$stage=="iv"]<-"III+IV"
table(fff$grade)

n<-length(names(df));n
for (i in 1:6) {
  #i=1
  a2<-fff[,c(i,n)]
  a2[a2==""]<-NA
  a2<-na.omit(a2)
  names(a2)
  name2<-names(a2)[1]
  names(a2)[1]<-"y1"
  y<-data.frame(compare_means(score~y1, data=a2))
  library(dplyr)
  y1<-y
  my_comparisons<-list()
  for (t in 1:length(rownames(y1))) {
    compare<-c(y1$group1[t],y1$group2[t])
    my_comparisons[[t]]<-compare
  }
  library(ggpubr)
  table(a2$y1)
  p1<-ggplot(a2, aes(x=a2[,1],y=a2[,2],fill=a2[,1])) +
    geom_boxplot(width=0.6,notch=T,notchwidth=0.3)+ #notchwidth越小则越往里凹
    labs(title=name2, x="", y="Risk score")+
    theme_bw()+theme(legend.position="none")+
    stat_compare_means(comparisons=my_comparisons,label="p.signif",
                       label.x = 1.5)
  pdf(file = paste("Signature_",name2,".pdf",sep = ""),width =2,height = 4,onefile = F);print(p1);dev.off()
}
