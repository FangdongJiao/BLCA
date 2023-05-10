myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");
setwd(myfilepath)

###################################LASSO
{
    rm(list = ls())
    load("cox.Rdata")
    cox_group<-cox_data[,1:3]
    df<-df[cox_group$GSM,]
    Lasdata3=as.matrix(df[,out1$Characteristics])#x自变量
    time2 = cbind(time = cox_group$time, status = cox_group$status)  #y应变量 
    set.seed(2022)
    library(glmnet)
    fit = glmnet(Lasdata3, time2, family="cox")#构建模型
    print(fit)
    dev.off()
    pdf("cox_lambda.pdf",width=8,height=5)
    par(las=0.5,mar=c(5,5,2,3))
    plot(fit, xvar="lambda", label=F)
    dev.off()
    
    cvfit = cv.glmnet(x=Lasdata3, y=time2, family = "cox")#进行交叉验证
    cvfit
    cvfit$lambda.min
    pdf("logistic_lambda.pdf",width=8,height=5)
    plot(cvfit)
    abline(v =log(cvfit$lambda.min) , col = "orange",lwd=2)
    dev.off()
    
    coefficients<-coef(cvfit, s ="lambda.min")#"lambda.min"
    Active.Index<-which(coefficients!=0)     #系数不为0的特征索引
    Active.coefficients<-coefficients[Active.Index]   #系数不为0的特征系数值
    gene<-coefficients@Dimnames[[1]][Active.Index]
    names(Active.coefficients)<-gene
    Active.coefficients1<-data.frame(Active.coefficients)
    save(Active.coefficients1,cox_group,df,file = "lasso_rick_score.Rdata")
  
}

##################################KM_ROC_风险因子联动图
{
  rm(list = ls())
  load("cox.Rdata")
  load("lasso_rick_score.Rdata")
  library(xlsx);write.xlsx(Active.coefficients1,"lasso_coef_exp.xlsx",
                           row.names = T)
  
  nams<-intersect(colnames(df),rownames(Active.coefficients1))
  expre_8<-df[,nams]
  data<-data.frame(expre_8)
  back_data<-cbind(cox_group,data)
  
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
  
  mgus1<-back_data
  library(timeROC)
  library(survival )
  ROC<-timeROC(T=mgus1$time,
               delta=mgus1$status,
               marker=mgus1$score,
               cause=1,
               weighting="marginal",
               times=c(365*1,365*3,365*5),
               ROC = TRUE,
               iid = TRUE)
  ROC
  
  library(survivalROC)
  require(ggsci)
  library("scales")
  dev.off()
  pdf("TCGA_ROC.pdf",width = 6,height = 6,onefile = F)
  plot(ROC,time=365*1,col = "red",add =FALSE, xlim=c(0,1), ylim=c(0,1),title="")#
  plot(ROC,time=365*3,col = "blue",add =T, xlim=c(0,1), ylim=c(0,1))#
  plot(ROC,time=365*5,col = "orange",add =T, xlim=c(0,1), ylim=c(0,1))#
  abline(0,1,col="black",lty=2)##线条颜色
  ROC_value<-data.frame(ROC$AUC)
  legend(0.4,0.2,c(paste("AUC of 1-year survival =",round(ROC_value[1,1],3)),
                   paste("AUC of 3-year survival =",round(ROC_value[2,1],3)),
                   paste("AUC of 5-year survival =",round(ROC_value[3,1],3))),
         x.intersp=1, y.intersp=0.8, #调整lengend在图中的位置
         lty= 1 ,lwd= 2,col=c("red","blue",'orange'),
         bty = "n",# bty框的类型
         seg.len=2,cex=0.8)
  dev.off()
  #KM
  fit_data<-back_data[,c(1:3,length(colnames(back_data)))]
  # 
  # ###################p>0.05  ,使用surv_cutpoint函数找到最优cutoff
  # res.cut <- surv_cutpoint(fit_data,
  #                          time = "time",
  #                          event = "status",
  #                          variables = "score",
  #                          minprop = 0.3) #可以添加多列
  # summary(res.cut)#查看最佳cutoff
  
  fit_data$levle<-ifelse(fit_data[,4]>median(fit_data[,4]),'high','low')
  table(fit_data$levle)
  KM_data<-fit_data
  library(survival)
  library(survminer)
  library(grid)
  fit<-survfit(Surv(time, status)~levle, data=KM_data)
  p<-  ggsurvplot(fit, data = KM_data, surv.median.line = "hv",
                  xlab="Time(days)",
                  palette =c("red","blue"),
                  legend.labs=c("high", "low"),
                  legend.title="TCGA_BLCA",
                  conf.int = F,risk.table = TRUE, pval=TRUE);p
  dev.off()
  pdf("TCGA_BLCA_KM_data.pdf",width=5,height=6,onefile = F);p;dev.off()
  
  library(xlsx);write.xlsx(KM_data,"TCGA_BLCA_risk_score.xlsx",row.names = F)
  save(KM_data,file = "TCGA_KM_data.Rdata")
  
  load("TCGA_KM_data.Rdata")
  #C风险因子联动图-表达量
  library(rms)
  library(ggrisk)
  RiskScore<-fit_data$score
  names(RiskScore) = rownames(fit_data)
  
  #开始绘制风险模型的生存点图
  fp <- RiskScore
  phe<-fit_data
  fp_dat=data.frame(patientid=1:length(fp),fp=as.numeric(sort(fp)))
  
  #添加风险分组，以风险评分的中位值将患者分为两组，大于中位值的 患者为高风险组，小于或等于中位值的患者为低风 险组
  fp_dat$riskgroup= ifelse(fp_dat$fp>median(fp_dat$fp),'high','low')
  
  sur_dat=data.frame(patientid=1:length(fp),
                     time=phe[names(fp),'time'],
                     event=phe[names(fp),'status'])
  sur_dat$event=ifelse(sur_dat$event==0,'alive','death')
  sur_dat$event=factor(sur_dat$event,levels = c("death","alive"))
  
  
  #fp_dat用来绘制第一幅图
  #sur_dat用来绘制第二幅图
  #exp_dat用来绘制第三幅图
  
  ###第一个图
  library(ggplot2)
  p1=ggplot(fp_dat,aes(x=patientid,y=fp))+geom_point(aes(color=riskgroup))+
    scale_colour_manual(values = c("red","BLUE"))+
    theme_bw()+labs(x="",y="Risk score")+
    geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8);p1
  dev.off()
  pdf("TCGA-风险因子联动图A.pdf",width=6,height=2,onefile=F);p1;dev.off()
  # png("TCGA-LIHC风险因子联动图A.png");p1;dev.off()
  #第二个图
  p2=ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=event))+theme_bw()+
    scale_colour_manual(values = c("red","BLUE"))+
    labs(x=" ",y="Survival time(days)")+
    geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8);p2
  dev.off()
  pdf("TCGA-风险因子联动图b.pdf",width=6,height=2,onefile=F);p2;  dev.off()
  # png("TCGA-LIHC风险因子联动图b.png"); p2; dev.off()
  
  # library(pheatmap)
  # mycolors <- colorRampPalette(c("BLUE", "WHITE", "red"), bias = 1.2)(100)
  # tmp=t(scale(data))
  # tmp[tmp > 1] = 1
  # tmp[tmp < -1] = -1
  # 
  # group_sample=KM_data[,c(1,5)]
  # group_sample<-group_sample[order(group_sample$levle),]
  # group_sample<-group_sample[,-1]
  # group_sample<-data.frame(group_sample)
  # rownames(group_sample)=KM_data$GSM
  # df<-group_sample
  # names(df)<-"group"
  # library(ComplexHeatmap)
  # table(df$group)
  # ha = HeatmapAnnotation(df=df,
  #                        col = list(group  = c("high"="red", 
  #                                              "low"="blue")))
  # p3<-Heatmap(tmp,cluster_rows = T, cluster_columns = F,
  #             col = colorRampPalette(c("red", "white", "blue")) (255),
  #             show_row_names = T,show_column_names = F,
  #             name = " ",
  #             top_annotation = ha);p3
  # dev.off()
  # # pdf("风险因子联动图c.pdf",width=4.5,height=3,onefile=F);p3;  dev.off()
  # #拼图实现三图联动
  # library(ggplotify)
  # plots = list(p1,p2,as.ggplot(as.grob(p3)))
  # library(gridExtra)
  # dev.off()
  # pdf("TCGA_风险因子联动图.pdf",width=10,height=15,onefile=F)
  # lay1 = rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7))) 
  # grid.arrange(grobs = plots, layout_matrix = lay1,
  #              heigths = c(2,3,2),weights=c(10,10,10))
  # dev.off()
}

###################################复杂热图
{
  rm(list = ls())
  load("lasso_rick_score.Rdata")
  heatmap_data<-df[,rownames(Active.coefficients1)]
  heatmap_data<-data.frame(t(heatmap_data))
  colnames(heatmap_data)<-gsub("\\.","-",colnames(heatmap_data))

  load("TCGA_KM_data.Rdata")
  table(KM_data$levle)
  KM_data<-na.omit(KM_data)
  KM_data<-KM_data[,c(1,5)]
  
  load("TCGA-BLCA.Rdata")
  data<-merge(KM_data,df1,by="GSM")
  df<-data[,-c(5,8)]
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
  fff[is.na(fff)]<-NA
 
  library(ComplexHeatmap)
  rownames(fff)<-fff$GSM
  data<-fff[,-1]
  names(data)
  ha  = HeatmapAnnotation(
    #"Group"=factor(data[,1]),
                          "Age"=factor(data[,2]),
                          "Gender"=factor(data[,3]),
                          "lymphatic_invasion"=factor(data[,4]),
                          "Stage"=factor(data[,5]),
                          "grade"=factor(data[,6]),
                          "therapy_outcome"=factor(data[,7]),
                          border = FALSE, annotation_name_side ="left",
                          simple_anno_size_adjust = TRUE
  )
  data<-data[order(data$levle,decreasing = F),]
  mat<-heatmap_data[,rownames(data)]
  
  tmp=t(scale(t(heatmap_data)))
  tmp[tmp > 1.5] = 1.5
  tmp[tmp < -1.5] = -1.5
  
  ha1 = HeatmapAnnotation(levle=data$levle,
                          col = list(levle  = c("high"="#fbb034", 
                                                  "low"="#0cb9c1")))
  hmap <-Heatmap(tmp,
                # column_km = 2,
                name = " ",
                cluster_columns = F,
                 show_row_names = T,
                 show_column_names = F,
                 # row_names_side = 'left',
                 bottom_annotation = ha,
                 top_annotation=ha1
  ); hmap
  dev.off()
  
  pdf("LASSO_heatmap.pdf",width = 10,height = 6,onefile = F)
  draw(hmap,
       heatmap_legend_side = 'left',
       annotation_legend_side = 'right')
  dev.off()
}
