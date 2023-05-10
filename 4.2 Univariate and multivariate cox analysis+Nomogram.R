myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");
setwd(myfilepath)

#TCGA-单+多
{
  rm(list = ls())
  load("TCGA-BLCA.Rdata")
  load("TCGA_KM_data.Rdata")
  pdata<-df1
  char_df<-merge(KM_data,pdata,by="GSM")
  fff<-char_df[,c(1:3,5,6,7,10,12)]
  
  table(fff$age)
  fff$age<-ifelse(fff$age>60,">60","<=60")
  table(fff$stage)
  fff$stage[fff$stage=="i"]<-"I+II"
  fff$stage[fff$stage=="ii"]<-"I+II"
  fff$stage[fff$stage=="iii"]<-"III+IV"
  fff$stage[fff$stage=="iv"]<-"III+IV"
  table(fff$grade)
  fff$grade[fff$grade==""]<-NA
  
  phenotype<-fff
  
  fff$stage<-ifelse(fff$stage=="III+IV",2,1)
  table(fff$gender)
  fff$gender<-ifelse(fff$gender=="female",2,1)
  table(fff$levle)
  fff$levle<-ifelse(fff$levle=="high",2,1)

  fff$grade<-ifelse(fff$grade=="High Grade",2,1)
  
  library(survival)
  library(forestplot)
  library(stringr)
  y<-colnames(fff)[-c(1:3)]
  y
  surv<-Surv(time = fff$time,event = fff$status)
  fff$surv<-with(fff,surv)
  d<-data.frame()
  for (i in 1:length(y)) {
    x<-y[i]
    FML<-as.formula(paste0("surv~",x))
    Cox<-coxph(FML,data = fff) 
    Sum<-summary(Cox)
    Pvalue<-round(Sum$coefficients[,5],8)
    HR<-round(Sum$coefficients[,2],2)
    low<-round(Sum$coefficients[,3],3)
    upper<-round(Sum$coefficients[,4],3)
    CI<-paste0(round(Sum$conf.int[,3:4],2),collapse = "-") 
    
    Unicox<-data.frame("Characteristics"=x,
                       "Hazard Ratio"=HR,
                       "CI95"=CI,
                       "Upper"=upper,
                       "Lower"=low,
                       "P value"=Pvalue)
    d<-rbind(d,Unicox)
  }
  d
  unicox<-d
  
  fff$score<-char_df$score
  
  d<-fff[,-c(9:10)]
  y<-colnames(d)[-(1:3)]
  y
  fml <- as.formula(paste0('Surv(time, status)~', paste(y, sep = '', collapse = '+')))
  fml
  res.cox <- coxph(fml, data = fff)
  x <- summary(res.cox)
  x$coefficients
  x$conf.int
  mCIL<-round(x$conf.int[,3],2)
  mCIU<-round(x$conf.int[,4],2)
  mCI<-paste0(mCIL, "-",mCIU) 
  mPvalue<-x$coefficients[,5]
  mHR<-round(x$coefficients[,2],2)
  D2<-data.frame("Characteristics"=paste("d",1:length(mHR)),
                 "Hazard Ratio"=mHR,
                 "CI95"=mCI,
                 "Upper"=mCIU,
                 "Lower"=mCIL,
                 "P value"=mPvalue)
  D2
  mulcox<-D2
  
  save(unicox,mulcox,fff,phenotype,file = "BLCA_cox.Rda")
  
  rm(list = ls())
  load("BLCA_cox.Rda")
  table(fff$stage)
  unicox
  mulcox
  unicox$d<-c("(high vs low)",
              "(>60 vs <=60)",
              "(High vs Grade)",
              "(III+IV vs I+II)",
              "(male vs female)")
  unicox
  result<-unicox[,c(1,7,6,3,2,4,5)] 
  result
  ins <- function(x) {c(x, rep(NA, ncol(result)-1))}
  result[,1]<-""
  
  head(result1)
  
  result1<-rbind(c("Characteristics","p.value","CI95",NA,NA,NA),
                 ins("Risk score"),
                 result[1,c(2:7)],
                 ins("Age"),
                 result[2, c(2:7)], 
                 ins("Grade"),
                 result[3, c(2:7)], 
                 ins("Stage"),
                 result[4, c(2:7)], 
                 ins("Gender"),
                 result[5, c(2:7)], 
                 c(NA,NA,NA,NA,NA,NA)
  )
  result1
  result1$t<-ifelse(result1$Upper=='NA',T,F)
  result1$t[is.na(result1$t)]<-T
  dev.off()
  library(forestplot)
  pdf("TCGA_unicox.pdf",width=6,height=6,onefile = F)
  forestplot(result1[,c(1:3)],
             mean=result1[,4],
             lower=result1[,5],
             upper=result1[,6], 
             zero=1,boxsize=0.2,graph.pos=4,
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=2),
                             "12"= gpar(lwd=2,lty=1)),
             graphwidth = unit(.25,"npc"),
             xlab=" ",xticks=c(0,1,1.5,2,3,4) ,
             is.summary=result1$t,
             txt_gp=fpTxtGp(label=gpar(cex=1),ticks=gpar(cex=1),
                            xlab=gpar(cex=1.5),title=gpar(cex=1)),
             lwd.zero=1,lwd.ci=1.5,lwd.xaxis=2,lty.ci=1,
             ci.vertices=TRUE,ci.vertices.height=0.2,
             clip=c(0.1,8),ineheight=unit(8, 'mm'), 
             line.margin=unit(8,'mm'),colgap=unit(6,'mm'),fn.ci_norm="fpDrawDiamondCI", 
             title="TCGA Unicox",
             col=fpColors(box ='#021eaa',lines ='#021eaa',zero = "black")
  )
  dev.off()
  
  mulcox$d<-c("(high vs low)",
              "(>60 vs <=60)",
              "(High vs Grade)",
              "(III+IV vs I+II)",
              "(male vs female)")
  mulcox
  result<-mulcox[,c(1,7,6,3,2,4,5)] 
  result
  ins <- function(x) {c(x, rep(NA, ncol(result)-1))}
  result[,1]<-""
  result1<-rbind(c("Characteristics","p.value","CI95",NA,NA,NA),
                 ins("Risk score"),
                 result[1,c(2:7)],
                 ins("Age"),
                 result[2, c(2:7)], 
                 ins("Grade"),
                 result[3, c(2:7)], 
                 ins("Stage"),
                 result[4, c(2:7)], 
                 ins("Gender"),
                 result[5, c(2:7)], 
                 c(NA,NA,NA,NA,NA,NA)
  )
  result1
  result1$t<-ifelse(result1$Upper=='NA',T,F)
  result1$t[is.na(result1$t)]<-T
  dev.off()
  library(forestplot)
  pdf("TCGA_Mulcox.pdf",width=7,height=7,onefile = F)
  forestplot(result1[,c(1:3)],
             mean=result1[,4],
             lower=result1[,5],
             upper=result1[,6], 
             zero=1,boxsize=0.2,graph.pos=4,
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=2),
                             "12"= gpar(lwd=2,lty=1)),
             graphwidth = unit(.25,"npc"),
             xlab=" ",xticks=c(0,1,1.5,2,3,4) ,
             is.summary=result1$t,
             txt_gp=fpTxtGp(label=gpar(cex=1),ticks=gpar(cex=1),
                            xlab=gpar(cex=1.5),title=gpar(cex=1)),
             lwd.zero=1,lwd.ci=1.5,lwd.xaxis=2,lty.ci=1,
             ci.vertices=TRUE,ci.vertices.height=0.2,
             clip=c(0.1,8),ineheight=unit(8, 'mm'), 
             line.margin=unit(8,'mm'),colgap=unit(6,'mm'),fn.ci_norm="fpDrawDiamondCI", 
             title="TCGA Mulcox",
             col=fpColors(box ='#021eaa',lines ='#021eaa',zero = "black")
  )
  dev.off()
  }

#TCGA--Nomogram
{
  rm(list = ls())
  load("BLCA_cox.Rda")
  data<-phenotype[,-4]#high low risk score
  data$score<-fff$score
  
  library(rms)
  dd=datadist(data)
  options(datadist="dd")
  library(Hmisc); 
  library(grid); 
  library(lattice);
  library(Formula); 
  library(ggplot2) 
  library(rms)
  library(survival)
  
  dd=datadist(data);options(datadist="dd") 
  rownames(data)<-data$GSM
  data<-data[,-1]
  #建立COX回归方程
  names(data)
  #何根据患者的年龄age，stage，grade等来预测患者的死亡情况
  f2 <- psm(Surv(time,status) ~ age+stage+grade+gender+score,data =  data, dist='lognormal') 
  surv <- Survival(f2) # 构建生存概率函数
  nom <- nomogram(f2, fun=list(function(x) surv(365*1, x),
                               function(x) surv(365*3, x),
                               function(x) surv(365*5, x)),
                  funlabel=c("1-year Survival Probability",
                             "3-year Survival Probability",
                             "5-year Survival Probability"))
  dev.off()

  pdf("列线图.pdf",width=15,height=10)
  plot(nom, xfrac=.6)
  dev.off()
  
  f1 <- psm(Surv(time, status) ~  age+stage+grade+gender+score, 
            x=T, y=T, 
            dist='lognormal',
            data=data,)
  cal1 <- calibrate(f1, cmethod="KM", 
                    method="boot", 
                    u=360*1, 
                    B=228)
  pdf("1-year.pdf",width = 6,height = 5,onefile = F)
  plot(cal1,
       xlim=c(0,1),
       ylim=c(0,1)
  )
  dev.off()
  
  cal2 <- calibrate(f1, cmethod="KM", method="boot", u=360*3,B=228)
  pdf("3-year.pdf",width = 6,height = 5,onefile = F)
  plot(cal2,
       xlim=c(0,1),
       ylim=c(0,1))
  dev.off()
  
  cal3 <- calibrate(f1, cmethod="KM", method="boot", u=360*5,B=228)
  pdf("5-year.pdf",width = 6,height = 5,onefile = F)
  plot(cal3,
       xlim=c(0,1),
       ylim=c(0,1))
  dev.off()

}

#GSE13507-单+多
{
  rm(list = ls())
  load("GSE13507.Rdata")
  load("GSE13507_KM_data.Rdata")
  names(GSE13507_pdata)
  df1<-data.frame(GSE13507_pdata$geo_accession,
                  GSE13507_pdata$`AGE:ch1`,
                  GSE13507_pdata$`stage:ch1`,
                  GSE13507_pdata$`grade:ch1`,
                  GSE13507_pdata$`SEX:ch1`)
  names(df1)<-c("GSM","age","stage","grade","gender")
  table(df1$age)
  table(df1$stage)
  for (i in 1:nrow(df1)) {
    df1$stage[i]<- paste(unlist(strsplit(df1$stage[i],""))[1:2],collapse = "")
  }
  df1$stage[df1$stage=="Ta"]<-NA
  table(df1$stage)
  
  pdata<-df1
  names(KM_data)[1]<-"GSM"
  char_df<-merge(KM_data,pdata,by="GSM")
  fff<-char_df[,-4]
  
  table(fff$age)
  fff$age<-ifelse(fff$age>60,">60","<=60")
  table(fff$stage)
  fff$stage[fff$stage=="T1"]<-"T1-2"
  fff$stage[fff$stage=="T2"]<-"T1-2"
  fff$stage[fff$stage=="T3"]<-"T3-4"
  fff$stage[fff$stage=="T4"]<-"T3-4"
  table(fff$grade)
  
  phenotype<-fff
  
  fff$stage<-ifelse(fff$stage=="T3-4",2,1)
  table(fff$gender)
  fff$gender<-ifelse(fff$gender=="F",2,1)
  table(fff$levle)
  fff$levle<-ifelse(fff$levle=="high",2,1)
  table(fff$grade)
  fff$grade<-ifelse(fff$grade=="high",2,1)
  
  library(survival)
  library(forestplot)
  library(stringr)
  y<-colnames(fff)[-c(1:3)]
  y
  surv<-Surv(time = fff$time,event = fff$status)
  fff$surv<-with(fff,surv)
  d<-data.frame()
  for (i in 1:length(y)) {
    x<-y[i]
    FML<-as.formula(paste0("surv~",x))
    Cox<-coxph(FML,data = fff) 
    Sum<-summary(Cox)
    Pvalue<-round(Sum$coefficients[,5],8)
    HR<-round(Sum$coefficients[,2],2)
    low<-round(Sum$coefficients[,3],3)
    upper<-round(Sum$coefficients[,4],3)
    CI<-paste0(round(Sum$conf.int[,3:4],2),collapse = "-") 
    
    Unicox<-data.frame("Characteristics"=x,
                       "Hazard Ratio"=HR,
                       "CI95"=CI,
                       "Upper"=upper,
                       "Lower"=low,
                       "P value"=Pvalue)
    d<-rbind(d,Unicox)
  }
  d
  unicox<-d
  
  fff$score<-char_df$score
  
  d<-fff[,-c(9:10)]
  y<-colnames(d)[-(1:3)]
  y
  fml <- as.formula(paste0('Surv(time, status)~', paste(y, sep = '', collapse = '+')))
  fml
  res.cox <- coxph(fml, data = fff)
  x <- summary(res.cox)
  x$coefficients
  x$conf.int
  mCIL<-round(x$conf.int[,3],2)
  mCIU<-round(x$conf.int[,4],2)
  mCI<-paste0(mCIL, "-",mCIU) 
  mPvalue<-x$coefficients[,5]
  mHR<-round(x$coefficients[,2],2)
  D2<-data.frame("Characteristics"=paste("d",1:length(mHR)),
                 "Hazard Ratio"=mHR,
                 "CI95"=mCI,
                 "Upper"=mCIU,
                 "Lower"=mCIL,
                 "P value"=mPvalue)
  D2
  mulcox<-D2
  
  save(unicox,mulcox,fff,phenotype,file = "GSE13507_cox.Rda")
  
  rm(list = ls())
  load("GSE13507_cox.Rda")
  table(fff$stage)
  unicox
  unicox$d<-c("(high vs low)",
              "(>60 vs <=60)",
              "(T3-4 vs T1-2)",
              "(High vs Grade)",
              "(male vs female)")
  unicox
  result<-unicox[,c(1,7,6,3,2,4,5)] 
  result
  ins <- function(x) {c(x, rep(NA, ncol(result)-1))}
  result[,1]<-""
  
  result1<-rbind(c("Characteristics","p.value","CI95",NA,NA,NA),
                 ins("Risk score"),
                 result[1,c(2:7)],
                 ins("Age"),
                 result[2, c(2:7)], 
                 ins("Stage"),
                 result[3, c(2:7)], 
                 ins("Grade"),
                 result[4, c(2:7)], 
                 ins("Gender"),
                 result[5, c(2:7)], 
                 c(NA,NA,NA,NA,NA,NA)
  )
  result1
  result1$t<-ifelse(result1$Upper=='NA',T,F)
  result1$t[is.na(result1$t)]<-T
  dev.off()
  library(forestplot)
  pdf("GSE13507_unicox.pdf",width=6,height=6,onefile = F)
  forestplot(result1[,c(1:3)],
             mean=result1[,4],
             lower=result1[,5],
             upper=result1[,6], 
             zero=1,boxsize=0.2,graph.pos=4,
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=2),
                             "12"= gpar(lwd=2,lty=1)),
             graphwidth = unit(.25,"npc"),
             xlab=" ",xticks=c(0,1,1.5,2,3,4) ,
             is.summary=result1$t,
             txt_gp=fpTxtGp(label=gpar(cex=1),ticks=gpar(cex=1),
                            xlab=gpar(cex=1.5),title=gpar(cex=1)),
             lwd.zero=1,lwd.ci=1.5,lwd.xaxis=2,lty.ci=1,
             ci.vertices=TRUE,ci.vertices.height=0.2,
             clip=c(0.1,8),ineheight=unit(8, 'mm'), 
             line.margin=unit(8,'mm'),colgap=unit(6,'mm'),fn.ci_norm="fpDrawDiamondCI", 
             title="GSE13507 Unicox",
             col=fpColors(box ='#021eaa',lines ='#021eaa',zero = "black")
  )
  dev.off()
  
  mulcox$d<-c("(high vs low)",
              "(>60 vs <=60)",
              "(T3-4 vs T1-2)",
              "(High vs Grade)",
              "(male vs female)")
  mulcox
  result<-mulcox[,c(1,7,6,3,2,4,5)] 
  result
  ins <- function(x) {c(x, rep(NA, ncol(result)-1))}
  result[,1]<-""
  result1<-rbind(c("Characteristics","p.value","CI95",NA,NA,NA),
                 ins("Risk score"),
                 result[1,c(2:7)],
                 ins("Age"),
                 result[2, c(2:7)], 
                 ins("Stage"),
                 result[3, c(2:7)], 
                 ins("Grade"),
                 result[4, c(2:7)], 
                 ins("Gender"),
                 result[5, c(2:7)], 
                 c(NA,NA,NA,NA,NA,NA)
  )
  result1
  result1$t<-ifelse(result1$Upper=='NA',T,F)
  result1$t[is.na(result1$t)]<-T
  dev.off()
  library(forestplot)
  pdf("GSE13507_Mulcox.pdf",width=7,height=7,onefile = F)
  forestplot(result1[,c(1:3)],
             mean=result1[,4],
             lower=result1[,5],
             upper=result1[,6], 
             zero=1,boxsize=0.2,graph.pos=4,
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=2),
                             "12"= gpar(lwd=2,lty=1)),
             graphwidth = unit(.25,"npc"),
             xlab=" ",xticks=c(0,1,1.5,2,3,4) ,
             is.summary=result1$t,
             txt_gp=fpTxtGp(label=gpar(cex=1),ticks=gpar(cex=1),
                            xlab=gpar(cex=1.5),title=gpar(cex=1)),
             lwd.zero=1,lwd.ci=1.5,lwd.xaxis=2,lty.ci=1,
             ci.vertices=TRUE,ci.vertices.height=0.2,
             clip=c(0.1,8),ineheight=unit(8, 'mm'), 
             line.margin=unit(8,'mm'),colgap=unit(6,'mm'),fn.ci_norm="fpDrawDiamondCI", 
             title="GSE13507 Mulcox",
             col=fpColors(box ='#021eaa',lines ='#021eaa',zero = "black")
  )
  dev.off()
}

#GSE48276-单+多
{
  rm(list = ls())
  load("GSE48276.Rdata")
  load("GSE48276_KM_data.Rdata")
  names(GSE48276_pdata)
  df1<-data.frame(GSE48276_pdata$geo_accession,
                  GSE48276_pdata$`age (at diagnosis):ch1`,
                  GSE48276_pdata$`pstage:ch1`,
                  # GSE48276_pdata$grade,
                  GSE48276_pdata$`gender:ch1`)
  names(df1)<-c("GSM","age","stage",
                #"grade",
                "gender")
  table(df1$age)
  table(df1$stage)
  for (i in 1:nrow(df1)) {
    df1$stage[i]<- paste(unlist(strsplit(df1$stage[i],""))[2:3],collapse = "")
  }
  df1$stage[df1$stage=="NANA"]<-NA
  table(df1$stage)
  
  pdata<-df1
  names(KM_data)[1]<-"GSM"
  char_df<-merge(KM_data,pdata,by="GSM")
  fff<-char_df[,-4]
  
  table(fff$age)
  fff$age<-ifelse(fff$age>60,">60","<=60")
  table(fff$stage)
  fff$stage[fff$stage=="T1"]<-"T1-2"
  fff$stage[fff$stage=="T0"]<-"T1-2"
  fff$stage[fff$stage=="T2"]<-"T1-2"
  fff$stage[fff$stage=="T3"]<-"T3-4"
  fff$stage[fff$stage=="T4"]<-"T3-4"
  table(fff$gender)
  
  phenotype<-fff
  
  names(fff)

  table(fff$age)
  fff$age<-ifelse(fff$age==">60",2,1)
  table(fff$stage)
  fff$stage<-ifelse(fff$stage=="T3-4",2,1)
  table(fff$gender)
  fff$gender<-ifelse(fff$gender=="F",2,1)
  table(fff$levle)
  fff$levle<-ifelse(fff$levle=="high",2,1)
  #table(fff$grade)
  #fff$grade<-ifelse(fff$grade=="high",2,1)
  
  library(survival)
  library(forestplot)
  library(stringr)
  y<-colnames(fff)[-c(1:3)]
  y
  surv<-Surv(time = fff$time,event = fff$status)
  fff$surv<-with(fff,surv)
  d<-data.frame()
  for (i in 1:length(y)) {
    x<-y[i]
    FML<-as.formula(paste0("surv~",x))
    Cox<-coxph(FML,data = fff) 
    Sum<-summary(Cox)
    Pvalue<-round(Sum$coefficients[,5],8)
    HR<-round(Sum$coefficients[,2],2)
    low<-round(Sum$coefficients[,3],3)
    upper<-round(Sum$coefficients[,4],3)
    CI<-paste0(round(Sum$conf.int[,3:4],2),collapse = "-") 
    
    Unicox<-data.frame("Characteristics"=x,
                       "Hazard Ratio"=HR,
                       "CI95"=CI,
                       "Upper"=upper,
                       "Lower"=low,
                       "P value"=Pvalue)
    d<-rbind(d,Unicox)
  }
  d
  unicox<-d
  
  
  fff$score<-char_df$score
  
  d<-fff[,-c(ncol(fff),ncol(fff)-1)]
  names(d)
  d<-d[,-c(5,7)]
  y<-colnames(d)[-(1:3)]
  y
  fml <- as.formula(paste0('Surv(time, status)~', paste(y, sep = '', collapse = '+')))
  fml
  res.cox <- coxph(fml, data = fff)
  x <- summary(res.cox)
  x$coefficients
  x$conf.int
  mCIL<-round(x$conf.int[,3],2)
  mCIU<-round(x$conf.int[,4],2)
  mCI<-paste0(mCIL, "-",mCIU) 
  mPvalue<-x$coefficients[,5]
  mHR<-round(x$coefficients[,2],2)
  D2<-data.frame("Characteristics"=paste("d",1:length(mHR)),
                 "Hazard Ratio"=mHR,
                 "CI95"=mCI,
                 "Upper"=mCIU,
                 "Lower"=mCIL,
                 "P value"=mPvalue)
  D2
  mulcox<-D2
  
  save(unicox,mulcox,fff,phenotype,file = "GSE48276_cox.Rda")
  
  rm(list = ls())
  load("GSE48276_cox.Rda")
  table(fff$stage)
  unicox
  unicox$d<-c("(high vs low)",
              "(>60 vs <=60)",
              "(T3-4 vs T1-2)",
              # "(High vs Grade)",
              "(male vs female)")
  unicox
  result<-unicox[,c(1,7,6,3,2,4,5)] 
  result
  ins <- function(x) {c(x, rep(NA, ncol(result)-1))}
  result[,1]<-""
  
  result1<-rbind(c("Characteristics","p.value","CI95",NA,NA,NA),
                 ins("Risk score"),
                 result[1,c(2:7)],
                 ins("Age"),
                 result[2, c(2:7)], 
                 ins("Stage"),
                 result[3, c(2:7)], 
                 # ins("Grade"),
                 # result[4, c(2:7)], 
                 ins("Gender"),
                 result[4, c(2:7)], 
                 c(NA,NA,NA,NA,NA,NA)
  )
  result1
  result1$t<-ifelse(result1$Upper=='NA',T,F)
  result1$t[is.na(result1$t)]<-T
  dev.off()
  library(forestplot)
  pdf("GSE48276_unicox.pdf",width=6,height=6,onefile = F)
  forestplot(result1[,c(1:3)],
             mean=result1[,4],
             lower=result1[,5],
             upper=result1[,6], 
             zero=1,boxsize=0.2,graph.pos=4,
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=2),
                             "10"= gpar(lwd=2,lty=1)),
             graphwidth = unit(.25,"npc"),
             xlab=" ",xticks=c(0,1,1.5,2,3,4) ,
             is.summary=result1$t,
             txt_gp=fpTxtGp(label=gpar(cex=1),ticks=gpar(cex=1),
                            xlab=gpar(cex=1.5),title=gpar(cex=1)),
             lwd.zero=1,lwd.ci=1.5,lwd.xaxis=2,lty.ci=1,
             ci.vertices=TRUE,ci.vertices.height=0.2,
             clip=c(0.1,8),ineheight=unit(8, 'mm'), 
             line.margin=unit(8,'mm'),colgap=unit(6,'mm'),fn.ci_norm="fpDrawDiamondCI", 
             title="GSE48276 Unicox",
             col=fpColors(box ='#021eaa',lines ='#021eaa',zero = "black")
  )
  dev.off()
  
  mulcox
  mulcox$d<-c("(high vs low)",
              # "(>60 vs <=60)",
              "(T3-4 vs T1-2)"
              # "(High vs Grade)",
              #"(male vs female)"
              )
  mulcox
  result<-mulcox[,c(1,7,6,3,2,4,5)] 
  result
  ins <- function(x) {c(x, rep(NA, ncol(result)-1))}
  result[,1]<-""
  result
  result1<-rbind(c("Characteristics","p.value","CI95",NA,NA,NA),
                 ins("Risk score"),
                 result[1,c(2:7)],
                 # ins("Age"),
                 # result[2, c(2:7)], 
                 ins("Stage"),
                 result[2, c(2:7)], 
                 # ins("Grade"),
                 # result[4, c(2:7)], 
                 # ins("Gender"),
                 # result[4, c(2:7)], 
                 c(NA,NA,NA,NA,NA,NA)
  )
  result1
  result1$t<-ifelse(result1$Upper=='NA',T,F)
  result1$t[is.na(result1$t)]<-T
  dev.off()
  
  library(forestplot)
  pdf("GSE48276_Mulcox.pdf",width=7,height=7,onefile = F)
  forestplot(result1[,c(1:3)],
             mean=result1[,4],
             lower=result1[,5],
             upper=result1[,6], 
             zero=1,boxsize=0.2,graph.pos=4,
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=2),
                             "6"= gpar(lwd=2,lty=1)),
             graphwidth = unit(.25,"npc"),
             xlab=" ",xticks=c(0,1,1.5,2,3,4) ,
             is.summary=result1$t,
             txt_gp=fpTxtGp(label=gpar(cex=1),ticks=gpar(cex=1),
                            xlab=gpar(cex=1.5),title=gpar(cex=1)),
             lwd.zero=1,lwd.ci=1.5,lwd.xaxis=2,lty.ci=1,
             ci.vertices=TRUE,ci.vertices.height=0.2,
             clip=c(0.1,8),ineheight=unit(8, 'mm'), 
             line.margin=unit(8,'mm'),colgap=unit(6,'mm'),fn.ci_norm="fpDrawDiamondCI", 
             title="GSE48276 Mulcox",
             col=fpColors(box ='#021eaa',lines ='#021eaa',zero = "black")
  )
  dev.off()
}

#GSE19915-单+多
{
  rm(list = ls())
  load("GSE19915.Rdata")
  load("GSE19915_KM_data.Rdata")
  names(GSE19915_pdata)
  df1<-data.frame(GSE19915_pdata$geo_accession,
                  # GSE19915_pdata$age,
                  GSE19915_pdata$`tumor stage:ch1`,
                  GSE19915_pdata$`tumor grade:ch1`)
  names(df1)<-c("GSM","stage",
                "grade")
  table(df1$stage)
  table(df1$grade)
  
  pdata<-df1
  names(KM_data)[1]<-"GSM"
  char_df<-merge(KM_data,pdata,by="GSM")
  fff<-char_df[,-5]
  
  table(fff$stage)
  fff$stage[fff$stage=="T1"]<-"T1-3"
  fff$stage[fff$stage=="T3"]<-"T1-3"
  table(fff$grade)
  fff$grade[fff$grade=="G1"]<-"G1-2"
  fff$grade[fff$grade=="G2"]<-"G1-2"
  
  phenotype<-fff
  
  table(fff$stage)
  fff$stage<-ifelse(fff$stage=="Ta",2,1)
  table(fff$grade)
  fff$grade<-ifelse(fff$grade=="G3",2,1)
  # table(fff$levle)
  # fff$levle<-ifelse(fff$levle=="high",2,1)
  
  library(survival)
  library(forestplot)
  library(stringr)
  y<-colnames(fff)[-c(1:3)]
  y
  surv<-Surv(time = fff$time,event = fff$status)
  fff$surv<-with(fff,surv)
  d<-data.frame()
  for (i in 1:length(y)) {
    x<-y[i]
    FML<-as.formula(paste0("surv~",x))
    Cox<-coxph(FML,data = fff) 
    Sum<-summary(Cox)
    Pvalue<-round(Sum$coefficients[,5],8)
    HR<-round(Sum$coefficients[,2],2)
    low<-round(Sum$coefficients[,3],3)
    upper<-round(Sum$coefficients[,4],3)
    CI<-paste0(round(Sum$conf.int[,3:4],2),collapse = "-") 
    
    Unicox<-data.frame("Characteristics"=x,
                       "Hazard Ratio"=HR,
                       "CI95"=CI,
                       "Upper"=upper,
                       "Lower"=low,
                       "P value"=Pvalue)
    d<-rbind(d,Unicox)
  }
  d
  unicox<-d
  
  fff$score<-char_df$score
  
  d<-fff[,-c(ncol(fff))]
  names(d)
  # d<-d[,-6]
  y<-colnames(d)[-(1:3)]
  y
  fml <- as.formula(paste0('Surv(time, status)~', paste(y, sep = '', collapse = '+')))
  fml
  res.cox <- coxph(fml, data = fff)
  x <- summary(res.cox)
  x$coefficients
  x$conf.int
  mCIL<-round(x$conf.int[,3],2)
  mCIU<-round(x$conf.int[,4],2)
  mCI<-paste0(mCIL, "-",mCIU) 
  mPvalue<-x$coefficients[,5]
  mHR<-round(x$coefficients[,2],2)
  D2<-data.frame("Characteristics"=paste("d",1:length(mHR)),
                 "Hazard Ratio"=mHR,
                 "CI95"=mCI,
                 "Upper"=mCIU,
                 "Lower"=mCIL,
                 "P value"=mPvalue)
  D2
  mulcox<-D2
  
  save(unicox,mulcox,fff,phenotype,file = "GSE19915_cox.Rda")
  
  rm(list = ls())
  load("GSE19915_cox.Rda")
  table(phenotype$grade)
  unicox
  unicox$d<-c("(high vs low)",
              "(Ta vs T1-3)",
              "(G3 vs G1-2)")
  unicox
  result<-unicox[,c(1,7,6,3,2,4,5)] 
  result
  ins <- function(x) {c(x, rep(NA, ncol(result)-1))}
  result[,1]<-""
  
  result1<-rbind(c("Characteristics","p.value","CI95",NA,NA,NA),
                 ins("Risk score"),
                 result[1,c(2:7)],
                 ins("Stage"),
                 result[2, c(2:7)], 
                 ins("Grade"),
                 result[3, c(2:7)], 
                 c(NA,NA,NA,NA,NA,NA)
  )
  result1
  result1$t<-ifelse(result1$Upper=='NA',T,F)
  result1$t[is.na(result1$t)]<-T
  dev.off()
  library(forestplot)
  pdf("GSE19915_unicox.pdf",width=6,height=6,onefile = F)
  forestplot(result1[,c(1:3)],
             mean=result1[,4],
             lower=result1[,5],
             upper=result1[,6], 
             zero=1,boxsize=0.2,graph.pos=4,
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=2),
                             "8"= gpar(lwd=2,lty=1)),
             graphwidth = unit(.25,"npc"),
             xlab=" ",xticks=c(0,1,5,8) ,
             is.summary=result1$t,
             txt_gp=fpTxtGp(label=gpar(cex=1),ticks=gpar(cex=1),
                            xlab=gpar(cex=1.5),title=gpar(cex=1)),
             lwd.zero=1,lwd.ci=1.5,lwd.xaxis=2,lty.ci=1,
             ci.vertices=TRUE,ci.vertices.height=0.2,
             clip=c(0.1,8),ineheight=unit(8, 'mm'), 
             line.margin=unit(8,'mm'),colgap=unit(6,'mm'),fn.ci_norm="fpDrawDiamondCI", 
             title="GSE19915 Unicox",
             col=fpColors(box ='#021eaa',lines ='#021eaa',zero = "black")
  )
  dev.off()
  
  mulcox$d<-c("(high vs low)",
              "(Ta vs T1-3)",
              "(G3 vs G1-2)")
  mulcox
  result<-mulcox[,c(1,7,6,3,2,4,5)] 
  result
  ins <- function(x) {c(x, rep(NA, ncol(result)-1))}
  result[,1]<-""
  result1<-rbind(c("Characteristics","p.value","CI95",NA,NA,NA),
                 ins("Risk score"),
                 result[1,c(2:7)],
                 ins("Stage"),
                 result[2, c(2:7)], 
                 ins("Grade"),
                 result[3, c(2:7)], 
                 c(NA,NA,NA,NA,NA,NA)
  )
  result1
  result1$t<-ifelse(result1$Upper=='NA',T,F)
  result1$t[is.na(result1$t)]<-T
  dev.off()
  library(forestplot)
  pdf("GSE19915_Mulcox.pdf",width=7,height=7,onefile = F)
  forestplot(result1[,c(1:3)],
             mean=result1[,4],
             lower=result1[,5],
             upper=result1[,6], 
             zero=1,boxsize=0.2,graph.pos=4,
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=2),
                             "8"= gpar(lwd=2,lty=1)),
             graphwidth = unit(.25,"npc"),
             xlab=" ",xticks=c(0,1,15,35) ,
             is.summary=result1$t,
             txt_gp=fpTxtGp(label=gpar(cex=1),ticks=gpar(cex=1),
                            xlab=gpar(cex=1.5),title=gpar(cex=1)),
             lwd.zero=1,lwd.ci=1.5,lwd.xaxis=2,lty.ci=1,
             ci.vertices=TRUE,ci.vertices.height=0.2,
             clip=c(0.1,8),ineheight=unit(8, 'mm'), 
             line.margin=unit(8,'mm'),colgap=unit(6,'mm'),fn.ci_norm="fpDrawDiamondCI", 
             title="GSE19915 Mulcox",
             col=fpColors(box ='#021eaa',lines ='#021eaa',zero = "black")
  )
  dev.off()
}