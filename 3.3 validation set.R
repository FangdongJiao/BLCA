myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");
setwd(myfilepath)

#GSE13507
{
  rm(list = ls()) 
  load("GSE13507.Rdata")
  names(GSE13507_pdata)
  cox_group<-GSE13507_pdata[,c(2,59,65,55,62,63,49)]
  head(cox_group)
  table(cox_group$`overall survival:ch1`)
  names(cox_group)<-c("GSM","status","time","grade","sex","stage","age")
  table(cox_group$status)
  cox_group$status<-ifelse(cox_group$status=="death",1,0)
  
  cox_group<-cox_group[!duplicated(cox_group$GSM),]
  cox_group<-na.omit(cox_group)
  
  rownames(cox_group)<-cox_group$GSM
  cox_group1<-cox_group[,1:3]
  
  load("lasso_rick_score.Rdata")
  nams<-intersect(rownames(GSE13507_matrix),rownames(Active.coefficients1));nams
  #nams<-nams[-2];nams
  expre_8<-GSE13507_matrix[nams,]
  expre_8<-t(expre_8)
  data<-data.frame(expre_8)
  
  f<-intersect(cox_group1$GSM,rownames(data))
  back<-cox_group1
  data$sample<-rownames(data)
  data<-data[f,]
  
  back<-back[f,]
  names(back)[1]<-"sample"
  back_data<-merge(back,data,by="sample")
  
  for (i in 1:length(nams)) {
    a<-names(back_data)[i+3]
    a1<-Active.coefficients1[a,]
    x<-a1
    back_data[,i+3]<-x*back_data[,i+3]
  }
  
  for (i in 1:length(rownames(back_data))) {
    x<-back_data[i,4:length(colnames(back_data))]
    x<-t(x)
    x<-as.numeric(x)
    x<-sum(x)
    back_data$score[i]<-x
  }
  names(back_data)[2:3]<-c("status","time")
  
  #C风险因子联动图-表达量
  fit_data<-back_data
  
  fit_data$time<-as.numeric(fit_data$time)
  fit_data$status<-as.numeric(fit_data$status)
  
  #中位值-----------------
  fit_data$levle<-ifelse(fit_data$score> median(fit_data$score),'high','low')
  table(fit_data$levle)
  n<-length(colnames(fit_data))
  KM_data<-fit_data[,c(1:3,n-1,n)];head(KM_data)
  library(survival)
  library(survminer)
  library(grid)
  fit<-survfit(Surv(time, status)~levle, data=KM_data)
  p<-ggsurvplot(fit, data = KM_data, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE13507",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  #找最佳阈值-------------
  res.cut <- surv_cutpoint(fit_data,
                           time = "time",
                           event = "status",
                           variables = "score",
                           minprop = 0.3)
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  fit <- survfit(Surv(time, status) ~score, data = res.cat)
  p<-ggsurvplot(fit, data = res.cat, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE13507",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  # #上下四分位----------------
  # km<-KM_data
  # down<-fivenum(km$score)[2]#下四分位
  # up<-fivenum(km$score)[4]#上四分位
  # x<-km[km$score<down,];x$levle<-"low"
  # y<-km[km$score>up,];y$levle<-"high"
  # km2<-rbind(x,y)
  # head(km2)
  # fit <- survfit(Surv(time, status) ~levle, data = km2)
  # p<-ggsurvplot(fit, data = km2, surv.median.line = "hv",
  #               xlab="Time(month)",legend.title="GSE13507",
  #               legend.labs=c("high", "low"),
  #               palette =c("red","blue"),
  #               conf.int = F,risk.table = TRUE, pval=TRUE);p
  # 
  # dev.off()
  # pdf("GSE13507_KM_data.pdf",width=5,height=6,onefile = F);p;  dev.off()
  # 
  table(res.cat$score)
  head(res.cat)
  KM_data$levle<-res.cat$score;head(KM_data)
  table(KM_data$level)
  library(xlsx);write.xlsx(KM_data[,-5],"GSE13507_risk_score.xlsx",row.names = F)
  save(KM_data,file = "GSE13507_KM_data.Rdata")
  
  load("GSE13507_KM_data.Rdata")
  mgus1<-KM_data
  library(timeROC) #选取score这个变量进行预后分析
  library(survival )
  mgus1$time<-as.numeric(mgus1$time)
  mgus1<-na.omit(mgus1)
  max(mgus1$time)
  ROC<-timeROC(T=mgus1$time,#结局时间
               delta=mgus1$status,#生存结局
               marker=mgus1$score,#预测变量
               # other_markers =mgus1$score,
               cause=1,#阳性结局赋值，比如死亡，复发的赋值
               weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
               times=c(12*1,12*3,12*5),#时间点，选取1-3年生存率,AUC<0.75，就3，5，7，10年的
               ROC = TRUE,
               iid = TRUE)
  ROC
  
  library(survivalROC)
  require(ggsci)
  library("scales")
  dev.off()
  pdf("GSE13507_ROC.pdf",width = 6,height = 6,onefile = F)
  plot(ROC,time=12*1,col = "red",add =FALSE, xlim=c(0,1), ylim=c(0,1),title="")#
  plot(ROC,time=12*3,col = "blue",add =T, xlim=c(0,1), ylim=c(0,1))#
  plot(ROC,time=12*5,col = "orange",add =T, xlim=c(0,1), ylim=c(0,1))#
  
  abline(0,1,col="black",lty=2)##线条颜色
  ROC_value<-data.frame(ROC$AUC)
  legend(0.4,0.2,c(paste("AUC of 1-year survival =",round(ROC_value[1,1],3)),
                   paste("AUC of 3-year survival =",round(ROC_value[2,1],3)),
                   paste("AUC of 5-year survival =",round(ROC_value[3,1],3))
                   
  ),
  x.intersp=1, y.intersp=0.8, #调整lengend在图中的位置
  lty= 1 ,lwd= 2,col=c("red","blue",'orange',"purple","pink"),
  bty = "n",# bty框的类型
  seg.len=2,cex=0.8)
  dev.off()
}

#GSE48276
{
  rm(list = ls()) 
  load("GSE48276.Rdata")
  names(GSE48276_pdata)
  cox_group<-GSE48276_pdata[,c(2,47,50,44,48,40)]
  head(cox_group)
  names(cox_group)<-c("GSM","status","time","sex","stage","age")
  table(cox_group$status)
  
  cox_group$status<-ifelse(cox_group$status=="censored",1,0)
  table(cox_group$status)
  
  cox_group<-cox_group[!duplicated(cox_group$GSM),]
  cox_group<-na.omit(cox_group)
  
  rownames(cox_group)<-cox_group$GSM
  cox_group1<-cox_group[,1:3]
  
  load("lasso_rick_score.Rdata")
  nams<-intersect(rownames(GSE48276_matrix),rownames(Active.coefficients1));nams
  #nams<-nams[-2];nams
  expre_8<-GSE48276_matrix[nams,]
  expre_8<-t(expre_8)
  data<-data.frame(expre_8)
  
  f<-intersect(cox_group1$GSM,rownames(data))
  back<-cox_group1
  data$sample<-rownames(data)
  data<-data[f,]
  
  back<-back[f,]
  names(back)[1]<-"sample"
  back_data<-merge(back,data,by="sample")
  
  for (i in 1:length(nams)) {
    a<-names(back_data)[i+3]
    a1<-Active.coefficients1[a,]
    x<-a1
    back_data[,i+3]<-x*back_data[,i+3]
  }
  
  for (i in 1:length(rownames(back_data))) {
    x<-back_data[i,4:length(colnames(back_data))]
    x<-t(x)
    x<-as.numeric(x)
    x<-sum(x)
    back_data$score[i]<-x
  }
  names(back_data)[2:3]<-c("status","time")
  
  #C风险因子联动图-表达量
  fit_data<-back_data
  
  fit_data$time<-as.numeric(fit_data$time)
  fit_data$status<-as.numeric(fit_data$status)
  
  fit_data$levle<-ifelse(fit_data$score> median(fit_data$score),'high','low')
  table(fit_data$levle)
  
  n<-length(colnames(fit_data))
  KM_data<-fit_data[,c(1:3,n-1,n)];head(KM_data)
  
  library(survival)
  library(survminer)
  library(grid)
  fit<-survfit(Surv(time, status)~levle, data=KM_data)
  p<-ggsurvplot(fit, data = KM_data, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE48276",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  #函数找最佳阈值
  res.cut <- surv_cutpoint(fit_data,
                           time = "time",
                           event = "status",
                           variables = "score",
                           minprop = 0.3)
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  fit <- survfit(Surv(time, status) ~score, data = res.cat)
  p<-ggsurvplot(fit, data = res.cat, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE48276",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  # #上下四分位----------------
  # km<-KM_data
  # down<-fivenum(km$score)[2]#下四分位
  # up<-fivenum(km$score)[4]#上四分位
  # x<-km[km$score<down,];x$levle<-"low"
  # y<-km[km$score>up,];y$levle<-"high"
  # km2<-rbind(x,y)
  # head(km2)
  # fit <- survfit(Surv(time, status) ~levle, data = km2)
  # p<-ggsurvplot(fit, data = km2, surv.median.line = "hv",
  #               xlab="Time(month)",legend.title="GSE48276",
  #               legend.labs=c("high", "low"),
  #               palette =c("red","blue"),
  #               conf.int = F,risk.table = TRUE, pval=TRUE);p
  # 
  # dev.off()
  dev.off()
  pdf("GSE48276_KM_data.pdf",width=5,height=6,onefile = F);p;  dev.off()
  
  table(res.cat$score)
  head(res.cat)
  KM_data$levle<-res.cat$score;head(KM_data)
  table(KM_data$levle)
  library(xlsx);write.xlsx(KM_data[,-5],"GSE48276_risk_score.xlsx",row.names = F)
  save(KM_data,file = "GSE48276_KM_data.Rdata")
  
  load("GSE48276_KM_data.Rdata")
  mgus1<-KM_data
  library(timeROC) #选取score这个变量进行预后分析
  library(survival )
  mgus1$time<-as.numeric(mgus1$time)
  mgus1<-na.omit(mgus1)
  max(mgus1$time)
  ROC<-timeROC(T=mgus1$time,#结局时间
               delta=mgus1$status,#生存结局
               marker=mgus1$score,#预测变量
               # other_markers =mgus1$score,
               cause=1,#阳性结局赋值，比如死亡，复发的赋值
               weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
               times=c(12*1,12*3,12*5,12*7,12*10),#时间点，选取1-3年生存率,AUC<0.75，就3，5，7，10年的
               ROC = TRUE,
               iid = TRUE)
  ROC
  
  library(survivalROC)
  require(ggsci)
  library("scales")
  dev.off()
  pdf("GSE48276_ROC.pdf",width = 6,height = 6,onefile = F)
  plot(ROC,time=12*3,col = "red",add =FALSE, xlim=c(0,1), ylim=c(0,1),title="")#
  plot(ROC,time=12*5,col = "blue",add =T, xlim=c(0,1), ylim=c(0,1))#
  plot(ROC,time=12*10,col = "orange",add =T, xlim=c(0,1), ylim=c(0,1))#
  
  abline(0,1,col="black",lty=2)##线条颜色
  ROC_value<-data.frame(ROC$AUC)
  legend(0.4,0.2,c(paste("AUC of 3-year survival =",round(ROC_value[2,1],3)),
                   paste("AUC of 5-year survival =",round(ROC_value[3,1],3)),
                   paste("AUC of 10-year survival =",round(ROC_value[5,1],3))
                   
  ),
  x.intersp=1, y.intersp=0.8, #调整lengend在图中的位置
  lty= 1 ,lwd= 2,col=c("red","blue",'orange',"purple","pink"),
  bty = "n",# bty框的类型
  seg.len=2,cex=0.8)
  dev.off()
}

#GSE19915✔
{
  rm(list = ls()) 
  load("GSE19915.Rdata")
  GSE19915_matrix<-na.omit(GSE19915_matrix)
  
  names(GSE19915_pdata)
  cox_group<-GSE19915_pdata[,c(2,49,51,
                              56,55)]
  head(cox_group)
  names(cox_group)
  names(cox_group)<-c("GSM","status","time","stage","grade")
  table(cox_group$status)

  cox_group$status<-ifelse(cox_group$status=="YES",1,0)
  
  cox_group<-na.omit(cox_group)
  
  rownames(cox_group)<-cox_group$GSM
  cox_group1<-cox_group[,1:3]
  
  load("lasso_rick_score.Rdata")
  nams<-intersect(rownames(GSE19915_matrix),rownames(Active.coefficients1));nams
  #nams<-nams[-2];nams
  expre_8<-GSE19915_matrix[nams,]
  expre_8<-t(expre_8)
  data<-data.frame(expre_8)
  
  f<-intersect(cox_group1$GSM,rownames(data))
  back<-cox_group1
  data$sample<-rownames(data)
  data<-data[f,]
  
  back<-back[f,]
  names(back)[1]<-"sample"
  back_data<-merge(back,data,by="sample")
  
  for (i in 1:length(nams)) {
    a<-names(back_data)[i+3]
    a1<-Active.coefficients1[a,]
    x<-a1
    back_data[,i+3]<-x*back_data[,i+3]
  }
  
  for (i in 1:length(rownames(back_data))) {
    x<-back_data[i,4:length(colnames(back_data))]
    x<-t(x)
    x<-as.numeric(x)
    x<-sum(x)
    back_data$score[i]<-x
  }
  names(back_data)[2:3]<-c("status","time")
  
  #C风险因子联动图-表达量
  fit_data<-back_data
  
  fit_data$time<-as.numeric(fit_data$time)
  fit_data$status<-as.numeric(fit_data$status)
  
  fit_data$levle<-ifelse(fit_data$score> median(fit_data$score),'high','low')
  table(fit_data$levle)
  
  n<-length(colnames(fit_data))
  KM_data<-fit_data[,c(1:3,n-1,n)];head(KM_data)
  
  library(survival)
  library(survminer)
  library(grid)
  fit<-survfit(Surv(time, status)~levle, data=KM_data)
  p<-ggsurvplot(fit, data = KM_data, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE19915",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  #函数找最佳阈值
  res.cut <- surv_cutpoint(fit_data,
                           time = "time",
                           event = "status",
                           variables = "score",
                           minprop = 0.3)
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  fit <- survfit(Surv(time, status) ~score, data = res.cat)
  p<-ggsurvplot(fit, data = res.cat, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE19915",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  dev.off()
  pdf("GSE19915_KM_data.pdf",width=5,height=6,onefile = F);p;  dev.off()
  
  table(res.cat$score)
  head(res.cat)
  KM_data$levle<-res.cat$score;head(KM_data)
  table(KM_data$level)
  library(xlsx);write.xlsx(KM_data[,-5],"GSE19915_risk_score.xlsx",row.names = F)
  save(KM_data,file = "GSE19915_KM_data.Rdata")
  
  load("GSE19915_KM_data.Rdata")
  mgus1<-KM_data
  library(timeROC) #选取score这个变量进行预后分析
  library(survival )
  mgus1$time<-as.numeric(mgus1$time)
  mgus1<-na.omit(mgus1)
  max(mgus1$time)
  ROC<-timeROC(T=mgus1$time,#结局时间
               delta=mgus1$status,#生存结局
               marker=mgus1$score,#预测变量
               # other_markers =mgus1$score,
               cause=1,#阳性结局赋值，比如死亡，复发的赋值
               weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
               times=c(12*1,12*3,12*5,12*7,12*10),#时间点，选取1-3年生存率,AUC<0.75，就3，5，7，10年的
               ROC = TRUE,
               iid = TRUE)
  ROC
  
  library(survivalROC)
  require(ggsci)
  library("scales")
  dev.off()
  pdf("GSE19915_ROC.pdf",width = 6,height = 6,onefile = F)
  plot(ROC,time=12*1,col = "red",add =FALSE, xlim=c(0,1), ylim=c(0,1),title="")#
  plot(ROC,time=12*3,col = "blue",add =T, xlim=c(0,1), ylim=c(0,1))#
  plot(ROC,time=12*5,col = "orange",add =T, xlim=c(0,1), ylim=c(0,1))#
  
  abline(0,1,col="black",lty=2)##线条颜色
  ROC_value<-data.frame(ROC$AUC)
  legend(0.4,0.2,c(paste("AUC of 1-year survival =",round(ROC_value[1,1],3)),
                   paste("AUC of 3-year survival =",round(ROC_value[2,1],3)),
                   paste("AUC of 5-year survival =",round(ROC_value[3,1],3))
                   
  ),
  x.intersp=1, y.intersp=0.8, #调整lengend在图中的位置
  lty= 1 ,lwd= 2,col=c("red","blue",'orange',"purple","pink"),
  bty = "n",# bty框的类型
  seg.len=2,cex=0.8)
  dev.off()
}

#GSE48075 ×
{
  rm(list = ls()) 
  load("GSE48075.Rdata")
  GSE48075_matrix<-na.omit(GSE48075_matrix)
  
  names(GSE48075_pdata)
  cox_group<-GSE48075_pdata[,c(2,43,47,
                               38,42)]
  head(cox_group)
  names(cox_group)
  names(cox_group)<-c("GSM","status","time","age","gender")
  table(cox_group$status)
  cox_group$status<-ifelse(cox_group$status=="censored",1,0)
  
  rownames(cox_group)<-cox_group$GSM
  cox_group1<-cox_group[,1:3]
  cox_group1<-na.omit(cox_group1)
  
  load("lasso_rick_score.Rdata")
  nams<-intersect(rownames(GSE48075_matrix),rownames(Active.coefficients1));nams
  #nams<-nams[-2];nams
  expre_8<-GSE48075_matrix[nams,]
  expre_8<-t(expre_8)
  data<-data.frame(expre_8)
  
  f<-intersect(cox_group1$GSM,rownames(data))
  back<-cox_group1
  data$sample<-rownames(data)
  data<-data[f,]
  
  back<-back[f,]
  names(back)[1]<-"sample"
  back_data<-merge(back,data,by="sample")
  
  for (i in 1:length(nams)) {
    a<-names(back_data)[i+3]
    a1<-Active.coefficients1[a,]
    x<-a1
    back_data[,i+3]<-x*back_data[,i+3]
  }
  
  for (i in 1:length(rownames(back_data))) {
    x<-back_data[i,4:length(colnames(back_data))]
    x<-t(x)
    x<-as.numeric(x)
    x<-sum(x)
    back_data$score[i]<-x
  }
  names(back_data)[2:3]<-c("status","time")
  
  #C风险因子联动图-表达量
  fit_data<-back_data
  
  fit_data$time<-as.numeric(fit_data$time)
  fit_data$status<-as.numeric(fit_data$status)
  
  fit_data$levle<-ifelse(fit_data$score> median(fit_data$score),'high','low')
  table(fit_data$levle)
  
  n<-length(colnames(fit_data))
  KM_data<-fit_data[,c(1:3,n-1,n)];head(KM_data)
  
  library(survival)
  library(survminer)
  library(grid)
  # #上下四分位----------------
  # km<-KM_data
  # down<-fivenum(km$score)[2]#下四分位
  # up<-fivenum(km$score)[4]#上四分位
  # x<-km[km$score<down,];x$levle<-"low"
  # y<-km[km$score>up,];y$levle<-"high"
  # km2<-rbind(x,y)
  # head(km2)
  # fit <- survfit(Surv(time, status) ~levle, data = km2)
  # p<-ggsurvplot(fit, data = km2, surv.median.line = "hv",
  #               xlab="Time(month)",legend.title="GSE48075",
  #               legend.labs=c("high", "low"),
  #               palette =c("red","blue"),
  #               conf.int = F,risk.table = TRUE, pval=TRUE);p
  # 
  
  fit<-survfit(Surv(time, status)~levle, data=KM_data)
  p<-ggsurvplot(fit, data = KM_data, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE48075",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  #函数找最佳阈值
  res.cut <- surv_cutpoint(fit_data,
                           time = "time",
                           event = "status",
                           variables = "score",
                           minprop = 0.3)
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  fit <- survfit(Surv(time, status) ~score, data = res.cat)
  p<-ggsurvplot(fit, data = res.cat, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE48075",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  dev.off()
  pdf("GSE48075_KM_data.pdf",width=5,height=6,onefile = F);p;  dev.off()
  
  table(res.cat$score)
  head(res.cat)
  KM_data$levle<-res.cat$score;head(KM_data)
  table(KM_data$level)
  library(xlsx);write.xlsx(KM_data[,-5],"GSE48075_risk_score.xlsx",row.names = F)
  save(KM_data,file = "GSE48075_KM_data.Rdata")
  
  load("GSE48075_KM_data.Rdata")
  mgus1<-KM_data
  library(timeROC) #选取score这个变量进行预后分析
  library(survival )
  mgus1$time<-as.numeric(mgus1$time)
  mgus1<-na.omit(mgus1)
  max(mgus1$time)
  ROC<-timeROC(T=mgus1$time,#结局时间
               delta=mgus1$status,#生存结局
               marker=mgus1$score,#预测变量
               # other_markers =mgus1$score,
               cause=1,#阳性结局赋值，比如死亡，复发的赋值
               weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
               times=c(12*1,12*3,12*5,12*7,12*10),#时间点，选取1-3年生存率,AUC<0.75，就3，5，7，10年的
               ROC = TRUE,
               iid = TRUE)
  ROC
  
  library(survivalROC)
  require(ggsci)
  library("scales")
  dev.off()
  pdf("GSE48075_ROC.pdf",width = 6,height = 6,onefile = F)
  plot(ROC,time=12*1,col = "red",add =FALSE, xlim=c(0,1), ylim=c(0,1),title="")#
  plot(ROC,time=12*3,col = "blue",add =T, xlim=c(0,1), ylim=c(0,1))#
  plot(ROC,time=12*5,col = "orange",add =T, xlim=c(0,1), ylim=c(0,1))#
  
  abline(0,1,col="black",lty=2)##线条颜色
  ROC_value<-data.frame(ROC$AUC)
  legend(0.4,0.2,c(paste("AUC of 1-year survival =",round(ROC_value[1,1],3)),
                   paste("AUC of 3-year survival =",round(ROC_value[2,1],3)),
                   paste("AUC of 5-year survival =",round(ROC_value[3,1],3))
                   
  ),
  x.intersp=1, y.intersp=0.8, #调整lengend在图中的位置
  lty= 1 ,lwd= 2,col=c("red","blue",'orange',"purple","pink"),
  bty = "n",# bty框的类型
  seg.len=2,cex=0.8)
  dev.off()
}

#GSE69795 ×
{
  rm(list = ls()) 
  load("GSE69795.Rdata")
  GSE69795_matrix<-na.omit(GSE69795_matrix)
  
  names(GSE69795_pdata)
  cox_group<-GSE69795_pdata[,c(2,49,56,
                               38,42)]
  head(cox_group)
  names(cox_group)
  names(cox_group)<-c("GSM","status","time","age","gender")
  table(cox_group$status)
  cox_group$status<-ifelse(cox_group$status=="Censored",1,0)
  
  rownames(cox_group)<-cox_group$GSM
  cox_group1<-cox_group[,1:3]
  cox_group1<-na.omit(cox_group1)
  
  load("lasso_rick_score.Rdata")
  nams<-intersect(rownames(GSE69795_matrix),rownames(Active.coefficients1));nams
  #nams<-nams[-2];nams
  expre_8<-GSE69795_matrix[nams,]
  expre_8<-t(expre_8)
  data<-data.frame(expre_8)
  
  f<-intersect(cox_group1$GSM,rownames(data))
  back<-cox_group1
  data$sample<-rownames(data)
  data<-data[f,]
  
  back<-back[f,]
  names(back)[1]<-"sample"
  back_data<-merge(back,data,by="sample")
  
  for (i in 1:length(nams)) {
    a<-names(back_data)[i+3]
    a1<-Active.coefficients1[a,]
    x<-a1
    back_data[,i+3]<-x*back_data[,i+3]
  }
  
  for (i in 1:length(rownames(back_data))) {
    x<-back_data[i,4:length(colnames(back_data))]
    x<-t(x)
    x<-as.numeric(x)
    x<-sum(x)
    back_data$score[i]<-x
  }
  names(back_data)[2:3]<-c("status","time")
  
  #C风险因子联动图-表达量
  fit_data<-back_data
  
  fit_data$time<-as.numeric(fit_data$time)
  fit_data$status<-as.numeric(fit_data$status)
  
  fit_data$levle<-ifelse(fit_data$score> median(fit_data$score),'high','low')
  table(fit_data$levle)
  
  n<-length(colnames(fit_data))
  KM_data<-fit_data[,c(1:3,n-1,n)];head(KM_data)
  
  library(survival)
  library(survminer)
  library(grid)
  # #上下四分位----------------
  # km<-KM_data
  # down<-fivenum(km$score)[2]#下四分位
  # up<-fivenum(km$score)[4]#上四分位
  # x<-km[km$score<down,];x$levle<-"low"
  # y<-km[km$score>up,];y$levle<-"high"
  # km2<-rbind(x,y)
  # head(km2)
  # fit <- survfit(Surv(time, status) ~levle, data = km2)
  # p<-ggsurvplot(fit, data = km2, surv.median.line = "hv",
  #               xlab="Time(month)",legend.title=" ",
  #               legend.labs=c("high", "low"),
  #               palette =c("red","blue"),
  #               conf.int = F,risk.table = TRUE, pval=TRUE);p
  # 
  
  fit<-survfit(Surv(time, status)~levle, data=KM_data)
  p<-ggsurvplot(fit, data = KM_data, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE69795",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  #函数找最佳阈值
  res.cut <- surv_cutpoint(fit_data,
                           time = "time",
                           event = "status",
                           variables = "score",
                           minprop = 0.3)
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  fit <- survfit(Surv(time, status) ~score, data = res.cat)
  p<-ggsurvplot(fit, data = res.cat, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE69795",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  dev.off()
  pdf("GSE69795_KM_data.pdf",width=5,height=6,onefile = F);p;  dev.off()
  
  table(res.cat$score)
  head(res.cat)
  KM_data$levle<-res.cat$score;head(KM_data)
  table(KM_data$level)
  library(xlsx);write.xlsx(KM_data[,-5],"GSE69795_risk_score.xlsx",row.names = F)
  save(KM_data,file = "GSE69795_KM_data.Rdata")
  
  load("GSE69795_KM_data.Rdata")
  mgus1<-KM_data
  library(timeROC) #选取score这个变量进行预后分析
  library(survival )
  mgus1$time<-as.numeric(mgus1$time)
  mgus1<-na.omit(mgus1)
  max(mgus1$time)
  ROC<-timeROC(T=mgus1$time,#结局时间
               delta=mgus1$status,#生存结局
               marker=mgus1$score,#预测变量
               # other_markers =mgus1$score,
               cause=1,#阳性结局赋值，比如死亡，复发的赋值
               weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
               times=c(12*1,12*3,12*5,12*7,12*10),#时间点，选取1-3年生存率,AUC<0.75，就3，5，7，10年的
               ROC = TRUE,
               iid = TRUE)
  ROC
  
  library(survivalROC)
  require(ggsci)
  library("scales")
  dev.off()
  pdf("GSE69795_ROC.pdf",width = 6,height = 6,onefile = F)
  plot(ROC,time=12*1,col = "red",add =FALSE, xlim=c(0,1), ylim=c(0,1),title="")#
  plot(ROC,time=12*3,col = "blue",add =T, xlim=c(0,1), ylim=c(0,1))#
  plot(ROC,time=12*5,col = "orange",add =T, xlim=c(0,1), ylim=c(0,1))#
  
  abline(0,1,col="black",lty=2)##线条颜色
  ROC_value<-data.frame(ROC$AUC)
  legend(0.4,0.2,c(paste("AUC of 1-year survival =",round(ROC_value[1,1],3)),
                   paste("AUC of 3-year survival =",round(ROC_value[2,1],3)),
                   paste("AUC of 5-year survival =",round(ROC_value[3,1],3))
                   
  ),
  x.intersp=1, y.intersp=0.8, #调整lengend在图中的位置
  lty= 1 ,lwd= 2,col=c("red","blue",'orange',"purple","pink"),
  bty = "n",# bty框的类型
  seg.len=2,cex=0.8)
  dev.off()
}

#GSE31684 
{
  rm(list = ls()) 
  load("GSE31684.Rdata")
  GSE31684_matrix<-na.omit(GSE31684_matrix)
  
  names(GSE31684_pdata)
  cox_group<-GSE31684_pdata[,c(2,75,74,
                               38,42)]
  head(cox_group)
  names(cox_group)
  names(cox_group)<-c("GSM","status","time","age","gender")
  table(cox_group$status)
  # cox_group$status<-ifelse(cox_group$status=="Censored",1,0)
  
  rownames(cox_group)<-cox_group$GSM
  cox_group1<-cox_group[,1:3]
  cox_group1<-na.omit(cox_group1)
  
  load("lasso_rick_score.Rdata")
  nams<-intersect(rownames(GSE31684_matrix),rownames(Active.coefficients1));nams
  #nams<-nams[-2];nams
  expre_8<-GSE31684_matrix[nams,]
  expre_8<-t(expre_8)
  data<-data.frame(expre_8)
  
  f<-intersect(cox_group1$GSM,rownames(data))
  back<-cox_group1
  data$sample<-rownames(data)
  data<-data[f,]
  
  back<-back[f,]
  names(back)[1]<-"sample"
  back_data<-merge(back,data,by="sample")
  
  for (i in 1:length(nams)) {
    a<-names(back_data)[i+3]
    a1<-Active.coefficients1[a,]
    x<-a1
    back_data[,i+3]<-x*back_data[,i+3]
  }
  
  for (i in 1:length(rownames(back_data))) {
    x<-back_data[i,4:length(colnames(back_data))]
    x<-t(x)
    x<-as.numeric(x)
    x<-sum(x)
    back_data$score[i]<-x
  }
  names(back_data)[2:3]<-c("status","time")
  
  #C风险因子联动图-表达量
  fit_data<-back_data
  
  fit_data$time<-as.numeric(fit_data$time)
  fit_data$status<-as.numeric(fit_data$status)
  
  fit_data$levle<-ifelse(fit_data$score> median(fit_data$score),'high','low')
  table(fit_data$levle)
  
  n<-length(colnames(fit_data))
  KM_data<-fit_data[,c(1:3,n-1,n)];head(KM_data)
  
  library(survival)
  library(survminer)
  library(grid)
  # #上下四分位----------------
  # km<-KM_data
  # down<-fivenum(km$score)[2]#下四分位
  # up<-fivenum(km$score)[4]#上四分位
  # x<-km[km$score<down,];x$levle<-"low"
  # y<-km[km$score>up,];y$levle<-"high"
  # km2<-rbind(x,y)
  # head(km2)
  # fit <- survfit(Surv(time, status) ~levle, data = km2)
  # p<-ggsurvplot(fit, data = km2, surv.median.line = "hv",
  #               xlab="Time(month)",legend.title=" ",
  #               legend.labs=c("high", "low"),
  #               palette =c("red","blue"),
  #               conf.int = F,risk.table = TRUE, pval=TRUE);p
  # 
  
  fit<-survfit(Surv(time, status)~levle, data=KM_data)
  p<-ggsurvplot(fit, data = KM_data, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE31684",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  #函数找最佳阈值
  res.cut <- surv_cutpoint(fit_data,
                           time = "time",
                           event = "status",
                           variables = "score",
                           minprop = 0.3)
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  fit <- survfit(Surv(time, status) ~score, data = res.cat)
  p<-ggsurvplot(fit, data = res.cat, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE31684",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  dev.off()
  pdf("GSE31684_KM_data.pdf",width=5,height=6,onefile = F);p;  dev.off()
  
  table(res.cat$score)
  head(res.cat)
  KM_data$levle<-res.cat$score;head(KM_data)
  table(KM_data$level)
  library(xlsx);write.xlsx(KM_data[,-5],"GSE31684_risk_score.xlsx",row.names = F)
  save(KM_data,file = "GSE31684_KM_data.Rdata")
  
  load("GSE31684_KM_data.Rdata")
  mgus1<-KM_data
  library(timeROC) #选取score这个变量进行预后分析
  library(survival )
  mgus1$time<-as.numeric(mgus1$time)
  mgus1<-na.omit(mgus1)
  max(mgus1$time)
  ROC<-timeROC(T=mgus1$time,#结局时间
               delta=mgus1$status,#生存结局
               marker=mgus1$score,#预测变量
               # other_markers =mgus1$score,
               cause=1,#阳性结局赋值，比如死亡，复发的赋值
               weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
               times=c(12*1,12*3,12*5,12*7,12*10),#时间点，选取1-3年生存率,AUC<0.75，就3，5，7，10年的
               ROC = TRUE,
               iid = TRUE)
  ROC
  
  library(survivalROC)
  require(ggsci)
  library("scales")
  dev.off()
  pdf("GSE31684_ROC.pdf",width = 6,height = 6,onefile = F)
  plot(ROC,time=12*1,col = "red",add =FALSE, xlim=c(0,1), ylim=c(0,1),title="")#
  plot(ROC,time=12*3,col = "blue",add =T, xlim=c(0,1), ylim=c(0,1))#
  plot(ROC,time=12*5,col = "orange",add =T, xlim=c(0,1), ylim=c(0,1))#
  
  abline(0,1,col="black",lty=2)##线条颜色
  ROC_value<-data.frame(ROC$AUC)
  legend(0.4,0.2,c(paste("AUC of 1-year survival =",round(ROC_value[1,1],3)),
                   paste("AUC of 3-year survival =",round(ROC_value[2,1],3)),
                   paste("AUC of 5-year survival =",round(ROC_value[3,1],3))
                   
  ),
  x.intersp=1, y.intersp=0.8, #调整lengend在图中的位置
  lty= 1 ,lwd= 2,col=c("red","blue",'orange',"purple","pink"),
  bty = "n",# bty框的类型
  seg.len=2,cex=0.8)
  dev.off()
}

#GSE69795 
{
  rm(list = ls()) 
  load("GSE69795.Rdata")
  GSE69795_matrix<-na.omit(GSE69795_matrix)
  
  names(GSE69795_pdata)
  cox_group<-GSE69795_pdata[,c(2,49,56,
                               38,42)]
  head(cox_group)
  names(cox_group)
  names(cox_group)<-c("GSM","status","time","age","gender")
  table(cox_group$status)
  cox_group$status<-ifelse(cox_group$status=="Censored",0,1)
  #发生终点事件记为“1”，删失(Censored)记为“0”
  table(cox_group$status)
  
  rownames(cox_group)<-cox_group$GSM
  cox_group1<-cox_group[,1:3]
  cox_group1<-na.omit(cox_group1)
  
  load("lasso_rick_score.Rdata")
  nams<-intersect(rownames(GSE69795_matrix),rownames(Active.coefficients1));nams
  #nams<-nams[-2];nams
  expre_8<-GSE69795_matrix[nams,]
  expre_8<-t(expre_8)
  data<-data.frame(expre_8)
  
  f<-intersect(cox_group1$GSM,rownames(data))
  back<-cox_group1
  data$sample<-rownames(data)
  data<-data[f,]
  
  back<-back[f,]
  names(back)[1]<-"sample"
  back_data<-merge(back,data,by="sample")
  
  for (i in 1:length(nams)) {
    a<-names(back_data)[i+3]
    a1<-Active.coefficients1[a,]
    x<-a1
    back_data[,i+3]<-x*back_data[,i+3]
  }
  
  for (i in 1:length(rownames(back_data))) {
    x<-back_data[i,4:length(colnames(back_data))]
    x<-t(x)
    x<-as.numeric(x)
    x<-sum(x)
    back_data$score[i]<-x
  }
  names(back_data)[2:3]<-c("status","time")
  
  #C风险因子联动图-表达量
  fit_data<-back_data
  
  fit_data$time<-as.numeric(fit_data$time)
  fit_data$status<-as.numeric(fit_data$status)
  
  fit_data$levle<-ifelse(fit_data$score> median(fit_data$score),'high','low')
  table(fit_data$levle)
  
  n<-length(colnames(fit_data))
  KM_data<-fit_data[,c(1:3,n-1,n)];head(KM_data)
  
  library(survival)
  library(survminer)
  library(grid)
  # #上下四分位----------------
  # km<-KM_data
  # down<-fivenum(km$score)[2]#下四分位
  # up<-fivenum(km$score)[4]#上四分位
  # x<-km[km$score<down,];x$levle<-"low"
  # y<-km[km$score>up,];y$levle<-"high"
  # km2<-rbind(x,y)
  # head(km2)
  # fit <- survfit(Surv(time, status) ~levle, data = km2)
  # p<-ggsurvplot(fit, data = km2, surv.median.line = "hv",
  #               xlab="Time(month)",legend.title=" ",
  #               legend.labs=c("high", "low"),
  #               palette =c("red","blue"),
  #               conf.int = F,risk.table = TRUE, pval=TRUE);p
  # 
  
  fit<-survfit(Surv(time, status)~levle, data=KM_data)
  p<-ggsurvplot(fit, data = KM_data, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE69795",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  #函数找最佳阈值
  res.cut <- surv_cutpoint(fit_data,
                           time = "time",
                           event = "status",
                           variables = "score",
                           minprop = 0.3)
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  fit <- survfit(Surv(time, status) ~score, data = res.cat)
  p<-ggsurvplot(fit, data = res.cat, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE69795",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  dev.off()
  pdf("GSE69795_KM_data.pdf",width=5,height=6,onefile = F);p;  dev.off()
  
  table(res.cat$score)
  head(res.cat)
  KM_data$levle<-res.cat$score;head(KM_data)
  table(KM_data$level)
  library(xlsx);write.xlsx(KM_data[,-5],"GSE69795_risk_score.xlsx",row.names = F)
  save(KM_data,file = "GSE69795_KM_data.Rdata")
  
  load("GSE69795_KM_data.Rdata")
  mgus1<-KM_data
  library(timeROC) #选取score这个变量进行预后分析
  library(survival )
  mgus1$time<-as.numeric(mgus1$time)
  mgus1<-na.omit(mgus1)
  max(mgus1$time)
  
  ROC<-timeROC(T=mgus1$time,#结局时间
               delta=mgus1$status,#生存结局
               marker=mgus1$score,#预测变量
               # other_markers =mgus1$score,
               cause=1,#阳性结局赋值，比如死亡，复发的赋值
               weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
               times=c(12*1,12*3,12*5),#时间点，选取1-3年生存率,AUC<0.75，就3，5，7，10年的
               ROC = TRUE,
               iid = TRUE);ROC
  
  library(survivalROC)
  require(ggsci)
  library("scales")
  dev.off()
  pdf("GSE69795_ROC.pdf",width = 6,height = 6,onefile = F)
  plot(ROC,time=12*1,col = "red",add =FALSE, xlim=c(0,1), ylim=c(0,1),title="")#
  plot(ROC,time=12*3,col = "blue",add =T, xlim=c(0,1), ylim=c(0,1))#
  plot(ROC,time=12*5,col = "orange",add =T, xlim=c(0,1), ylim=c(0,1))#
  
  abline(0,1,col="black",lty=2)##线条颜色
  ROC_value<-data.frame(ROC$AUC)
  legend(0.4,0.2,c(paste("AUC of 1-year survival =",round(ROC_value[1,1],3)),
                   paste("AUC of 3-year survival =",round(ROC_value[2,1],3)),
                   paste("AUC of 5-year survival =",round(ROC_value[3,1],3))
                   
  ),
  x.intersp=1, y.intersp=0.8, #调整lengend在图中的位置
  lty= 1 ,lwd= 2,col=c("red","blue",'orange',"purple","pink"),
  bty = "n",# bty框的类型
  seg.len=2,cex=0.8)
  dev.off()
}

#GSE39280 
{
  rm(list = ls()) 
  load("GSE39280.Rdata")
  GSE39280_matrix<-na.omit(GSE39280_matrix)
  
  names(GSE39280_pdata)
  cox_group<-GSE39280_pdata[,c(2,83,64,
                               38,42)]
  head(cox_group)
  names(cox_group)
  names(cox_group)<-c("GSM","status","time","age","gender")
  table(cox_group$status)
  cox_group$status<-ifelse(cox_group$status=="living",0,1)
  #发生终点事件记为“1”，删失(Censored)记为“0”
  table(cox_group$status)
  
  rownames(cox_group)<-cox_group$GSM
  cox_group1<-cox_group[,1:3]
  cox_group1<-na.omit(cox_group1)
  
  load("lasso_rick_score.Rdata")
  nams<-intersect(rownames(GSE39280_matrix),rownames(Active.coefficients1));nams
  #nams<-nams[-2];nams
  expre_8<-GSE39280_matrix[nams,]
  expre_8<-t(expre_8)
  data<-data.frame(expre_8)
  
  f<-intersect(cox_group1$GSM,rownames(data))
  back<-cox_group1
  data$sample<-rownames(data)
  data<-data[f,]
  
  back<-back[f,]
  names(back)[1]<-"sample"
  back_data<-merge(back,data,by="sample")
  
  for (i in 1:length(nams)) {
    a<-names(back_data)[i+3]
    a1<-Active.coefficients1[a,]
    x<-a1
    back_data[,i+3]<-x*back_data[,i+3]
  }
  
  for (i in 1:length(rownames(back_data))) {
    x<-back_data[i,4:length(colnames(back_data))]
    x<-t(x)
    x<-as.numeric(x)
    x<-sum(x)
    back_data$score[i]<-x
  }
  names(back_data)[2:3]<-c("status","time")
  
  #C风险因子联动图-表达量
  fit_data<-back_data
  
  fit_data$time<-as.numeric(fit_data$time)
  fit_data$status<-as.numeric(fit_data$status)
  
  fit_data$levle<-ifelse(fit_data$score> median(fit_data$score),'high','low')
  table(fit_data$levle)
  
  n<-length(colnames(fit_data))
  KM_data<-fit_data[,c(1:3,n-1,n)];head(KM_data)
  
  library(survival)
  library(survminer)
  library(grid)
  # #上下四分位----------------
  # km<-KM_data
  # down<-fivenum(km$score)[2]#下四分位
  # up<-fivenum(km$score)[4]#上四分位
  # x<-km[km$score<down,];x$levle<-"low"
  # y<-km[km$score>up,];y$levle<-"high"
  # km2<-rbind(x,y)
  # head(km2)
  # fit <- survfit(Surv(time, status) ~levle, data = km2)
  # p<-ggsurvplot(fit, data = km2, surv.median.line = "hv",
  #               xlab="Time(month)",legend.title=" ",
  #               legend.labs=c("high", "low"),
  #               palette =c("red","blue"),
  #               conf.int = F,risk.table = TRUE, pval=TRUE);p
  # 
  
  fit<-survfit(Surv(time, status)~levle, data=KM_data)
  p<-ggsurvplot(fit, data = KM_data, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE39280",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  #函数找最佳阈值
  res.cut <- surv_cutpoint(fit_data,
                           time = "time",
                           event = "status",
                           variables = "score",
                           minprop = 0.3)
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  fit <- survfit(Surv(time, status) ~score, data = res.cat)
  p<-ggsurvplot(fit, data = res.cat, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE39280",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  dev.off()
  pdf("GSE39280_KM_data.pdf",width=5,height=6,onefile = F);p;  dev.off()
  
  table(res.cat$score)
  head(res.cat)
  KM_data$levle<-res.cat$score;head(KM_data)
  table(KM_data$level)
  library(xlsx);write.xlsx(KM_data[,-5],"GSE39280_risk_score.xlsx",row.names = F)
  save(KM_data,file = "GSE39280_KM_data.Rdata")
  
  load("GSE39280_KM_data.Rdata")
  mgus1<-KM_data
  library(timeROC) #选取score这个变量进行预后分析
  library(survival )
  mgus1$time<-as.numeric(mgus1$time)
  mgus1<-na.omit(mgus1)
  max(mgus1$time)
  
  ROC<-timeROC(T=mgus1$time,#结局时间
               delta=mgus1$status,#生存结局
               marker=mgus1$score,#预测变量
               # other_markers =mgus1$score,
               cause=1,#阳性结局赋值，比如死亡，复发的赋值
               weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
               times=c(12*1,12*3,12*5),#时间点，选取1-3年生存率,AUC<0.75，就3，5，7，10年的
               ROC = TRUE,
               iid = TRUE);ROC
  
  library(survivalROC)
  require(ggsci)
  library("scales")
  dev.off()
  pdf("GSE39280_ROC.pdf",width = 6,height = 6,onefile = F)
  plot(ROC,time=12*1,col = "red",add =FALSE, xlim=c(0,1), ylim=c(0,1),title="")#
  plot(ROC,time=12*3,col = "blue",add =T, xlim=c(0,1), ylim=c(0,1))#
  plot(ROC,time=12*5,col = "orange",add =T, xlim=c(0,1), ylim=c(0,1))#
  
  abline(0,1,col="black",lty=2)##线条颜色
  ROC_value<-data.frame(ROC$AUC)
  legend(0.4,0.2,c(paste("AUC of 1-year survival =",round(ROC_value[1,1],3)),
                   paste("AUC of 3-year survival =",round(ROC_value[2,1],3)),
                   paste("AUC of 5-year survival =",round(ROC_value[3,1],3))
                   
  ),
  x.intersp=1, y.intersp=0.8, #调整lengend在图中的位置
  lty= 1 ,lwd= 2,col=c("red","blue",'orange',"purple","pink"),
  bty = "n",# bty框的类型
  seg.len=2,cex=0.8)
  dev.off()
}

#GSE37816 
{
  rm(list = ls()) 
  load("GSE37816.Rdata")
  GSE37816_matrix<-na.omit(GSE37816_matrix)
  
  names(GSE37816_pdata)
  cox_group<-GSE37816_pdata[,c(2,36,37)]
  #  ,38,42)]
  head(cox_group)
  names(cox_group)
  names(cox_group)<-c("GSM","status","time")
  #,"age","gender")
  table(cox_group$status)
  cox_group$status<-ifelse(cox_group$status=="survival, 2:death): 1",1,0)
  #发生终点事件记为“1”，删失(Censored)记为“0”
  table(cox_group$status)
  
  rownames(cox_group)<-cox_group$GSM
  cox_group1<-cox_group[,1:3]
  cox_group1<-na.omit(cox_group1)
  
  load("lasso_rick_score.Rdata")
  nams<-intersect(rownames(GSE37816_matrix),rownames(Active.coefficients1));nams
  #nams<-nams[-2];nams
  expre_8<-GSE37816_matrix[nams,]
  expre_8<-t(expre_8)
  data<-data.frame(expre_8)
  
  f<-intersect(cox_group1$GSM,rownames(data))
  back<-cox_group1
  data$sample<-rownames(data)
  data<-data[f,]
  
  back<-back[f,]
  names(back)[1]<-"sample"
  back_data<-merge(back,data,by="sample")
  
  for (i in 1:length(nams)) {
    a<-names(back_data)[i+3]
    a1<-Active.coefficients1[a,]
    x<-a1
    back_data[,i+3]<-x*back_data[,i+3]
  }
  
  for (i in 1:length(rownames(back_data))) {
    x<-back_data[i,4:length(colnames(back_data))]
    x<-t(x)
    x<-as.numeric(x)
    x<-sum(x)
    back_data$score[i]<-x
  }
  names(back_data)[2:3]<-c("status","time")
  
  #C风险因子联动图-表达量
  fit_data<-back_data
  
  fit_data$time<-as.numeric(fit_data$time)
  fit_data$status<-as.numeric(fit_data$status)
  
  fit_data$levle<-ifelse(fit_data$score> median(fit_data$score),'high','low')
  table(fit_data$levle)
  
  n<-length(colnames(fit_data))
  KM_data<-fit_data[,c(1:3,n-1,n)];head(KM_data)
  
  library(survival)
  library(survminer)
  library(grid)
  # #上下四分位----------------
  # km<-KM_data
  # down<-fivenum(km$score)[2]#下四分位
  # up<-fivenum(km$score)[4]#上四分位
  # x<-km[km$score<down,];x$levle<-"low"
  # y<-km[km$score>up,];y$levle<-"high"
  # km2<-rbind(x,y)
  # head(km2)
  # fit <- survfit(Surv(time, status) ~levle, data = km2)
  # p<-ggsurvplot(fit, data = km2, surv.median.line = "hv",
  #               xlab="Time(month)",legend.title=" ",
  #               legend.labs=c("high", "low"),
  #               palette =c("red","blue"),
  #               conf.int = F,risk.table = TRUE, pval=TRUE);p
  # 
  
  fit<-survfit(Surv(time, status)~levle, data=KM_data)
  p<-ggsurvplot(fit, data = KM_data, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE37816",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  #函数找最佳阈值
  res.cut <- surv_cutpoint(fit_data,
                           time = "time",
                           event = "status",
                           variables = "score",
                           minprop = 0.3)
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  fit <- survfit(Surv(time, status) ~score, data = res.cat)
  p<-ggsurvplot(fit, data = res.cat, surv.median.line = "hv",
                xlab="Time(month)",legend.title="GSE37816",
                legend.labs=c("high", "low"),
                palette =c("red","blue"),
                conf.int = F,risk.table = TRUE, pval=TRUE);p
  
  dev.off()
  pdf("GSE37816_KM_data.pdf",width=5,height=6,onefile = F);p;  dev.off()
  
  table(res.cat$score)
  head(res.cat)
  KM_data$levle<-res.cat$score;head(KM_data)
  table(KM_data$level)
  library(xlsx);write.xlsx(KM_data[,-5],"GSE37816_risk_score.xlsx",row.names = F)
  save(KM_data,file = "GSE37816_KM_data.Rdata")
  
  load("GSE37816_KM_data.Rdata")
  mgus1<-KM_data
  library(timeROC) #选取score这个变量进行预后分析
  library(survival )
  mgus1$time<-as.numeric(mgus1$time)
  mgus1<-na.omit(mgus1)
  max(mgus1$time)
  
  ROC<-timeROC(T=mgus1$time,#结局时间
               delta=mgus1$status,#生存结局
               marker=mgus1$score,#预测变量
               # other_markers =mgus1$score,
               cause=1,#阳性结局赋值，比如死亡，复发的赋值
               weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
               times=c(12*1,12*3,12*5),#时间点，选取1-3年生存率,AUC<0.75，就3，5，7，10年的
               ROC = TRUE,
               iid = TRUE);ROC
  
  library(survivalROC)
  require(ggsci)
  library("scales")
  dev.off()
  pdf("GSE37816_ROC.pdf",width = 6,height = 6,onefile = F)
  plot(ROC,time=12*1,col = "red",add =FALSE, xlim=c(0,1), ylim=c(0,1),title="")#
  plot(ROC,time=12*3,col = "blue",add =T, xlim=c(0,1), ylim=c(0,1))#
  plot(ROC,time=12*5,col = "orange",add =T, xlim=c(0,1), ylim=c(0,1))#
  
  abline(0,1,col="black",lty=2)##线条颜色
  ROC_value<-data.frame(ROC$AUC)
  legend(0.4,0.2,c(paste("AUC of 1-year survival =",round(ROC_value[1,1],3)),
                   paste("AUC of 3-year survival =",round(ROC_value[2,1],3)),
                   paste("AUC of 5-year survival =",round(ROC_value[3,1],3))
                   
  ),
  x.intersp=1, y.intersp=0.8, #调整lengend在图中的位置
  lty= 1 ,lwd= 2,col=c("red","blue",'orange',"purple","pink"),
  bty = "n",# bty框的类型
  seg.len=2,cex=0.8)
  dev.off()
}