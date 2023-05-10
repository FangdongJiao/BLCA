rm(list = ls())
myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");myfilepath
setwd(myfilepath)

rm(list = ls())
library(oncoPredict)
dir='./IC50_DataFiles/Training Data'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
dim(GDSC2_Expr)  
GDSC2_Expr[1:4, 1:4]
GDSC2_Expr<-data.frame(GDSC2_Expr)

load("lasso_rick_score.Rdata")
#系数*样本表达值
d<-intersect(rownames(Active.coefficients1),rownames(GDSC2_Expr))
expre_8<-GDSC2_Expr[d,]
expre_8<-na.omit(expre_8)
expre_8<-t(expre_8)
data<-data.frame(expre_8)
back_data<-data
for (i in 1:length(colnames(back_data))) {
  # i=1
  a<-names(back_data)[i]
  a1<-Active.coefficients1[a,]
  x<-a1
  back_data[,i]<-x*back_data[,i]
}
#计算每组风险评分，（相加）
for (i in 1:length(rownames(back_data))) {
  x<-back_data[i,]
  x<-t(x)
  x<-as.numeric(x)
  x<-sum(x)
  back_data$score[i]<-x
}
back_data$GSM<-rownames(back_data)
fit_data<-back_data[,c(length(colnames(back_data)),length(colnames(back_data))-1)]
#中值判断高低风险组
fit_data$levle<-ifelse(fit_data[,2]>median(fit_data[,2]),'high','low')

# Read GDSC2 response data. 
#行——样本
#列——药物
dir
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res<-data.frame(GDSC2_Res)
dim(GDSC2_Res)  # 805 198
GDSC2_Res[1:4,1:4]
names(GDSC2_Res)=gsub('_[0-9]*','',colnames(GDSC2_Res))
GDSC2_df<-data.frame(t(GDSC2_Res))
rownames(GDSC2_df)<-tolower(rownames(GDSC2_df))

# x<-rownames(GDSC2_df)[1:50]
x<-c("cisplatin","gemcitabine","gefitinib","doxorubicin",
     "cediranib","carmustine","camptothecin","buparlisib",
     "afatinib")
df<-GDSC2_df[x,]
df<-df[-4,]
df1<-data.frame(t(df))
df1$GSM<-rownames(df1)
data<-merge(df1,fit_data,by="GSM")

library(reshape2)
data1<-melt(data)
#data1<-na.omit(data1)
data1[1:4,1:4]
table(unique(data1$variable))
data1<-data1[data1$variable!="score",]
table(unique(data1$variable))

library(ggpubr)
library(ggthemes)
library(ggbeeswarm)
p<-ggplot(data1,aes(x=levle,y=value))+geom_boxplot()+
  geom_beeswarm(aes(color=factor(levle)),cex=2,alpha = 0.3)+
  theme_bw()+ xlab("")+ylab("Esitimated IC50")+
  scale_color_manual(values=c("red","blue"))+guides(color=F)+
  stat_compare_means(aes(group =  levle),
                     label = "p.signif",label.x = 1.5,
                     hide.ns = F)+
  facet_wrap(~variable,scales="free");p
dev.off()  
pdf("Drug_signature.pdf",onefile = F,width = 15,height = 12);p;dev.off()  

#gene_IC50_hign-low------------------------------------
gene_IG50<-data.frame(GDSC2_Expr[d,])

dim(gene_IG50)  # 805 198
gene_expre<-data.frame(t(gene_IG50))
gene_expre$GSM<-rownames(gene_expre)
gh<-merge(gene_expre,fit_data,by="GSM")
gh1<-melt(gh)
head(gh1)
p<-ggplot(gh1,aes(x=levle,y=value))+geom_boxplot(aes(colour=levle))+
  theme_bw()+ xlab("")+ylab("Esitimated IC50")+
  scale_color_manual(values=c("red","blue"))+
  guides(color=F)+
  stat_compare_means(aes(group =  levle),
                     label = "p.signif",label.x = 1.5,
                     hide.ns = F)+
  facet_wrap(~variable,scales="free");p
dev.off()  
pdf("Gene_Drug_signature.pdf",onefile = F,width = 30,height = 24);p;dev.off()  

#df[is.na(df)]<-0
combination <- expand.grid(rownames(gene_IG50),rownames(df))
data<-rbind(gene_IG50,df)
names(combination)=c("Gene","Drug")
library(corrplot)
cor_result=apply(combination,1,function(x){
  drug=as.character(x[1])
  KEGG=as.character(x[2])
  result=cor.test(as.numeric(data[drug,]), as.numeric(data[KEGG,]),
                  method ="spearman")
  score=c(pval=result$p.value,result$estimate)
  return(score)
})
result=cbind(combination,t(cor_result))

head(result)
library(ggplot2)
library(stringr)
p1<-ggplot(result,aes(Drug,Gene))+
  geom_point(aes(size=pval,color=rho))+       
  scale_color_gradient2(low = "#0cb9c1", mid = "white", high = "#e4002b")+
  labs(color=expression(rho),size="-log10(pval)", x="",y="",title="")+ 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1));p1
dev.off()
pdf("Gene_drug_cor.pdf",width = 5,height = 8,onefile = F);p1;dev.off()
