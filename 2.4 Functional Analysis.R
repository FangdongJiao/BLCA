myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");
setwd(myfilepath)

rm(list = ls())
load("gene_modelu.Rda")

file<-list.files(pattern = "Hs.gmt")
geneset<-data.frame()
for (i in 1:length(file)) {
  # i=1
  o<- clusterProfiler::read.gmt(file[i]) 
  geneset<-rbind(geneset,o)
}
gene_set<-geneset[,c(2,1)]

names(GENE)
load("gene_modelu.Rda")
names(GENE)
B<-c(GENE[[1]]$gene,
     GENE[[2]]$gene
     # GENE[[3]]$gene
)
x<-intersect(B,geneset$gene)

library(org.Hs.eg.db);library(clusterProfiler)
library(ggplot2)
library(stringr)

ENTREZID <-mapIds(org.Hs.eg.db, 
                  keys=B,
                  column="ENTREZID", 
                  keytype="SYMBOL" )
ENTREZID<-unlist(ENTREZID)
ENTREZID<-data.frame(ENTREZID)
ENTREZID<-na.omit(ENTREZID)

go=enrichGO(ENTREZID$ENTREZID,
            OrgDb = org.Hs.eg.db,
            ont = 'ALL',
            pAdjustMethod = 'BH',
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.2,
            keyType = 'ENTREZID')
go<-filter(go,p.adjust<0.05)
dim(go)
go<-data.frame(go)
table(go$ONTOLOGY)

#install.packages('R.utils')
R.utils::setOption( "clusterProfiler.download.method",'auto' )
kegg <- enrichKEGG(ENTREZID$ENTREZID, 
                   organism = 'hsa', 
                   keyType = 'kegg',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH', 
                   minGSSize = 0,
                   maxGSSize = 100000000,
                   qvalueCutoff = 0.2,
                   use_internal_data = FALSE)
kegg<-filter(kegg,p.adjust<0.05)
dim(kegg)
KEGG<-data.frame(kegg)
save(go,KEGG,file = "Module_GO_KEGG.Rda")

rm(list = ls())
load("Module_GO_KEGG.Rda")
library(DOSE)
redio <- data.frame(parse_ratio(KEGG$GeneRatio))
Round<-round(redio$parse_ratio.KEGG.GeneRatio., 3)
Round<-data.frame(Round)
names(Round)<-"GeneRatio"
KEGG$GeneRatio<-Round$GeneRatio
KEGG<-KEGG[order(KEGG$Count,decreasing = T),]
head(KEGG)

p1 <- ggplot(KEGG[1:20,],
             aes(GeneRatio,
                 factor(Description,
                        levels =Description)))+
  geom_point(aes(size=GeneRatio,
                 color=p.adjust))+  
  scale_color_gradient(low = "blue", high = "red")+
  labs(color=expression(p.adjust),
       size="GeneRatio", 
       x="GeneRatio",
       y=" ",
       title="KEGG")+ 
  theme_bw() +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 60))+ #scale_y(x)_discrete,以防文本过长
  theme(plot.title = element_text(hjust =0.5));p1
dev.off()
pdf("Module_kegg.pdf",width = 8,height = 8,onefile = F);p1;dev.off()

library(DOSE)
redio <- data.frame(parse_ratio(go$GeneRatio))
Round<-round(redio$parse_ratio.go.GeneRatio., 3)
Round<-data.frame(Round)
names(Round)<-"GeneRatio"
go$GeneRatio<-Round$GeneRatio

BP=go[go$ONTOLOGY=='BP',]
CC=go[go$ONTOLOGY=='CC',]
MF=go[go$ONTOLOGY=='MF',]

go_data<-list("BP"=BP,"CC"=CC,"MF"=MF)
for (i in 1:length(go_data)) {
  # i=1
  df<-go_data[[i]]
  alldata<-df[order(df$Count,decreasing = T),];head(alldata)
  library(ggplot2)
  p1<-ggplot(alldata[1:20,],
             aes(GeneRatio,
                 factor(Description,
                        levels =Description)))+
    geom_point(aes(size=GeneRatio,
                   color=p.adjust))+  
    scale_color_gradient(low = "blue", high = "red")+
    labs(color=expression(p.adjust),
         size="GeneRatio", 
         x="GeneRatio",y="",
         title= names(go_data[i]))+ 
    theme_bw() +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 60))+ #scale_y(x)_discrete,以防文本过长
    theme(plot.title = element_text(hjust =0.5))
  pdf(paste("Module_go_",names(go_data[i]),".pdf",sep = ""),
      width = 8,height = 8,onefile = F);print(p1);dev.off()
}

