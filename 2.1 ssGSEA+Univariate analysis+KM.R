rm(list = ls())
myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");myfilepath
setwd(myfilepath)

#ssGSEA----------------------------------------------------------------------
rm(list = ls())
load("TCGA-BLCA.Rdata")
table(TCGA_group$group_list)

file<-list.files(pattern = "Hs.gmt")
geneset<-data.frame()
for (i in 1:length(file)) {
        # i=1
        o<- clusterProfiler::read.gmt(file[i]) 
        geneset<-rbind(geneset,o)
}
gene_set<-geneset[,c(2,1)]
list<- split(as.matrix(gene_set)[,1], gene_set[,2])

library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
gsva_matrix<- gsva(as.matrix(TCGA_exp), 
                   list,method='ssgsea',
                   kcdf='Gaussian',abs.ranking=TRUE)
save(gsva_matrix,TCGA_group,file = 'ssGSEA.Rdata')


#单因素--------------------------------------------------------------------------
rm(list = ls())
load("TCGA-BLCA.Rdata")
load("ssGSEA.Rdata")
matrix<-data.frame(t(gsva_matrix),check.rows = T)

library(survival)
matrix$sample<-rownames(matrix)

cox_data<-merge(survival,matrix,by="sample")

cox_data$OS<-as.numeric(cox_data$OS)
cox_data$OS.time<-as.numeric(cox_data$OS.time)

library(survminer)
library(tableone)  
library(survival)
memory.limit(size=100000)##60498  7863
y<-colnames(cox_data)[-c(1:3)]

resulit_all=data.frame()
for (i in 1:length(y)){
        candidate_gene<-y[i]
        data<-cox_data[,c(2,3,i+3)]
        formula <- as.formula(paste0('Surv(OS.time, OS)~', candidate_gene))
        surv_uni_cox <- coxph(formula, data = data)
        m<-data.frame(ShowRegTable(surv_uni_cox))#https://www.jianshu.com/p/52232599fc3b
        resulit_all=rbind (resulit_all,m)
}
out1<-resulit_all
out<-out1
library(stringr)
for (g in 1:length(out$exp.coef...confint.)) {
        out$HR[g]<-strsplit(out$exp.coef...confint.,"\\[")[[g]][1]
}
out$CI<-str_extract(out$exp.coef...confint.,"\\[.+?\\]")#提取中括号内的字符串
x<-gsub("\\[","",out$CI)
x1<-gsub("\\]","",x)
x2<-strsplit(x1,",")
for (i in 1:length(x2)) {
        out$lower[i]<-x2[[i]][1]
}
for (i in 1:length(x2)) {
        out$upper[i]<-x2[[i]][2]
}
out<-out[,-c(1,4)]
out<-out[,c(2:4,1)]
names(out)[2:3]<-c("lower95_CI","upper95_CI")
out$hallmaker<-rownames(out)
hallmaker_50_CI<-out[,c(5,1:4)]

save(cox_data,hallmaker_50_CI,file = "ssGSEA.cox.Rdata")

rm(list = ls())
load("ssGSEA.cox.Rdata")
sig_dif1<-hallmaker_50_CI
rownames <- cbind(c("Pathway" , sig_dif1$hallmaker),
                  c("Hazard_Ratio", sig_dif1$HR),
                  c("lower 95% CI", sig_dif1$lower95_CI),
                  c("upper 95% CI", sig_dif1$upper95_CI),
                  c("pvalue", sig_dif1$p))
library(dplyr)
df <- tibble(sig_dif1[,1:5])
test_data <- rbind(rep(NA, 5), df)
names(test_data)
sig_dif1<-rbind(c("Pathway", "HR","low 95%CI","upper 95%CI","p.value"),
                sig_dif1[1:4, ],
                c(NA,NA,NA,NA,NA,NA))
library(forestplot)
library(tibble)
pdf("Hallmark Univariate.pdf",width=15,height=5,onefile=F)
forestplot(sig_dif1,  #https://mp.weixin.qq.com/s/4lkSp8Qh3s8RDaS0ldqh1Q
           test_data,
           graph.pos=6,#p.value设置图出现在第几列
           mean=as.numeric(sig_dif1$HR),
           lower=as.numeric(sig_dif1$lower95_CI), 
           upper=as.numeric(sig_dif1$upper95_CI),
           boxsize=0.3,#方块大小
           zero=1, #x轴垂直线
           hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                           "2" = gpar(lty=2),
                           "6"= gpar(lwd=2,lty=1)),
           graphwidth = unit(.25,"npc"),
           xlab=" ",xticks=c(0,0.5,1) ,
           # is.summary=result1$t,
           txt_gp=fpTxtGp(label=gpar(cex=1),ticks=gpar(cex=1),
                          xlab=gpar(cex=1.5),title=gpar(cex=1)),
           lwd.zero=1,lwd.ci=1.5,lwd.xaxis=2,lty.ci=1,
           #箱线图两端添加小竖线，高度
           ci.vertices=TRUE,ci.vertices.height=0.2,
           clip=c(0.1,8),ineheight=unit(8, 'mm'), 
           line.margin=unit(8,'mm'),colgap=unit(6,'mm'),fn.ci_norm="fpDrawDiamondCI", 
           title=" ",
           col=fpColors(box ='#021eaa',lines ='#021eaa',zero = "black")
)
dev.off()

#KM---------------------------------------------------------------
rm(list = ls())
load("cox.Rdata")
load("TCGA-BLCA.Rdata")
table(TCGA_group$group_list)
tumor<-TCGA_group[TCGA_group$group_list=="tumor",]
cox_data<-cox_data[cox_data$sample%in%tumor$GSM,]

survival<-read.table("TCGA_survival.txt",sep = "\t",header = T)
survival$sample<-paste(survival$sample,"A",sep = "")
names(survival)
d<-cox_data[,-c(2,3)]

library(survival)
library(survminer)
library(grid)
y<-c(3,5,7,9) 
for (x in 1:4) {
        # x=1
        p<-y[x]        
        survival1<-survival[,c(1,p,p+1)]
        back_data<-merge(survival1,d,by="sample")
        for (i in 1:length(1:4)) {
                # i=1
                fit_data<-back_data[,c(1:3,i+3)]
                unit<-names(fit_data)[2:3]
                names(fit_data)[2:3]<-c("OS","OS.time")
                
                name<-names(fit_data)[4]
                
                ########################median
                fit_data$levle<-ifelse(fit_data[,4]>median(fit_data[,4]),'high','low')
                table(fit_data$levle)
                KM_data<-na.omit(fit_data)
               
                head(KM_data)
                fit<-survfit(Surv(OS.time, OS)~levle, data=KM_data)
                p<-  ggsurvplot(fit, data = KM_data, surv.median.line = "hv",
                                xlab=unit[1],
                                palette =c("red","blue"),
                                legend.labs=c("high", "low"),
                                legend.title=name,
                                conf.int = F,risk.table = F, pval=TRUE)
                pdf(paste("./2.1、KM曲线/",name,"_KM_",unit[1],".pdf",sep = ""),width=5,height=5,onefile = F);print(p);dev.off()
                
                #----------------------------
                fit_data<-back_data[,c(1:3,i+3)]
                unit<-names(fit_data)[2:3]
                names(fit_data)[2:3]<-c("OS","OS.time")
                
                name<-names(fit_data)[4];names(fit_data)[4]<-"score"
                res.cut <- surv_cutpoint(fit_data,
                                         time = "OS.time",
                                         event = "OS",
                                         variables = "score",
                                         minprop = 0.3)
                res.cat <- surv_categorize(res.cut)
                fit <- survfit(Surv(OS.time, OS) ~score, data = res.cat)
                p<-  ggsurvplot(fit, data = res.cat, surv.median.line = "hv",
                                xlab=unit[1],
                                palette =c("red","blue"),
                                legend.labs=c("high", "low"),
                                legend.title=name,
                                conf.int = F,risk.table = F, pval=TRUE);p
                pdf(paste("./2.1、KM曲线/",name,"_KM_",unit[1],"_surv_cutpoint.pdf",sep = ""),width=5,height=5,onefile = F);print(p);dev.off()
                
        }  
}


#boxplot_ssGSEA_score_high_low-------------------------------------------
rm(list = ls())
load("ssGSEA.cox.Rdata")
ssGSEA<-cox_data[,-c(2:3)]

load("TCGA-BLCA.Rdata")

df<-df1
table(df$age)
df$age<-ifelse(df$age>50,">50","<=50")
df$therapy_outcome[df$therapy_outcome=="Complete Remission/Response"]<-"CR"
df$therapy_outcome[df$therapy_outcome=="Partial Remission/Response"]<-"PR"
df$therapy_outcome[df$therapy_outcome=="Progressive Disease"]<-"PD"
df$therapy_outcome[df$therapy_outcome=="Stable Disease"]<-"SD"

fff<-df
fff$stage[fff$stage=="i"]<-"I"
fff$stage[fff$stage=="ii"]<-"II"
fff$stage[fff$stage=="iii"]<-"III"
fff$stage[fff$stage=="iv"]<-"IV"

df1<-na.omit(fff)
names(df1)
df1<-df1[,-c(4,7)]
dir.create("ssGSEA_pdata")
for (u in 2:5) {
        #u=2
        data<-ssGSEA[,c(1,u)]
        names(data)[1]<-"GSM"
        
        for (i in 1:6) {
                #i=1
                d<-df1[,c(1,i+1)]
                a2<-merge(data,d,by="GSM")
                a2[a2==""]<-NA
                a2<-na.omit(a2)
                
                q0<-names(a2)[2]
                names(a2)[2]<-"pathway"
                
                q1<-names(a2)[3]
                names(a2)[3]<-"y1"
                
                y<-data.frame(compare_means(pathway~y1, data=a2))
                library(dplyr)
                y1<-y
                my_comparisons<-list()
                for (t in 1:length(rownames(y1))) {
                        compare<-c(y1$group1[t],y1$group2[t])
                        my_comparisons[[t]]<-compare
                }
                library(ggbeeswarm)
                library(ggpubr)
                head(a2)
                p1<-ggplot(a2, aes(x=y1,y=pathway,fill=y1)) +
                        geom_boxplot()+ 
                        labs(title=q1, x="", y="ssGSEA score")+
                        theme_bw()+theme(legend.position="none")+
                        stat_compare_means(comparisons=my_comparisons,
                                           label="p.signif",
                                           label.x = 1.5);
                pdf(file = paste("ssGSEA_pdata/",names(data)[2],
                                 names(d)[2],".pdf",sep = "_"),
                    width = 4,height = 4,onefile = F)
                print(p1)
                dev.off()  
        }
}
