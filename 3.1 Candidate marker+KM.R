myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");
setwd(myfilepath)

#单因素cox------------------------------------------------------------------
rm(list = ls())
load("TCGA-BLCA.Rdata")
boxplot(TCGA_exp[1:50,1:20])

load("gene_modelu.Rda")
names(GENE)
B<-c(GENE[[1]]$gene,
    GENE[[2]]$gene
     # GENE[[3]]$gene
     )
length(B)
data<-TCGA_exp[B,]
names(survival)<-c("GSM","status","time")
nam<-intersect(survival$GSM,colnames(data))

df<-data.frame(t(data[,nam]))
df$GSM<-rownames(df)
cox_group<-survival[nam,]
cox_data<-merge(cox_group,df,by="GSM")#生存数据及表达矩阵的组合
cox_data$time<-as.numeric(cox_data$time)
cox_data$status<-as.numeric(cox_data$status)

library('survival')
library(tableone)
memory.limit(size=100000)##60498  7863
resulit_all=data.frame()
y<-colnames(cox_data)[-c(1:3)];y

gene_Unicox<-data.frame()

for (i in 1:length(y)){
  # i=131
  x<-y[i]
  candidate_gene<-y[i];candidate_gene
  w<-cox_data[,c(2,3,i+3)]
  formula <- as.formula(paste0('Surv(time, status)~', candidate_gene))
  surv_uni_cox <- coxph(formula, data = w)
  Sum<-summary(surv_uni_cox)
  Pvalue<-round(Sum$coefficients[,5],4)
  HR<-round(Sum$coefficients[,2],2)
  low<-round(Sum$conf.int[,3],3)
  upper<-round(Sum$conf.int[,4],3)
  CI<-paste0(round(Sum$conf.int[,3:4],2),collapse = "-")
  Unicox<-data.frame("Characteristics"=x,
                     "Hazard Ratio"=HR,
                     "CI95"=CI,
                     "Upper"=upper,
                     "Lower"=low,
                     "P value"=Pvalue)
  gene_Unicox<-rbind(gene_Unicox,Unicox)
}
out<-gene_Unicox
out1<-out[out$P.value<0.001,]
save(data,cox_group,cox_data,out,out1,df,file = "cox.Rdata")

#
rm(list = ls())
load("cox.Rdata")
load("lasso_rick_score.Rdata")
rownames(out1)<-out1$Characteristics
gg<-out1[rownames(Active.coefficients1),]
f<-df[,rownames(gg)]
matrix<-data.frame(f)
matrix$GSM<-rownames(matrix)
cox_group<-cox_data[,1:3]
KM<-merge(cox_group,matrix,by="GSM")#生存数据及表达矩阵的组合
KM$time<-as.numeric(KM$time)

library(survival)
library(survminer) 
dir.create("LASSO_KM")
for (i in 4:length(names(KM))) {
  #i=7
  KM_data<-KM[,c(1:3,i)]
  KM_data$expression<-ifelse(KM_data[,4]>median(KM_data[,4]),'>median','<=median')##中位数判断大小
  fit<-survfit(Surv(time, status)~expression, data=KM_data)
  p3<-ggsurvplot(fit,  pval=F, #https://www.51xxziyuan.com/58/1202.html
                 conf.int = F,
                 # break.x.by = 1000,
                 xlab="Time(days)", palette = c("red", "blue"),##更改线条颜色
                 legend.labs=c("<=median", ">median"),
                 legend.title =names(KM_data[4]),
                 ggtheme = theme_bw())
  #https://www.136.la/jingpin/show-50803.html
  ###添加HR ,CI ,P
  res_cox<-coxph(Surv(time, status) ~expression, data=KM_data)
  p3$plot = p3$plot + 
    ggplot2::annotate("text",x = 1500, y = 0.1,
                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + 
    ggplot2::annotate("text",x = 1500, y = 0.15,
                      label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
    ggplot2::annotate("text",x = 1500, y = 0.2,
                      #size = 1.5,
                      label = paste("P:",round(summary(res_cox)$coef[5],5)))
  pdf(paste("LASSO_KM/",names(KM_data)[4],"_KM_",".pdf"), width = 4, height = 4,onefile = F)
  print(p3)
  dev.off()
  
  
}
