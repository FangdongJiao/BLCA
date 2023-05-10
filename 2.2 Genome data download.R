rm(list = ls())
myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");myfilepath
setwd(myfilepath)

rm(list = ls()) 
library(TCGAbiolinks)
query <- GDCquery(
  project = "TCGA-BLCA", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)#https://mp.weixin.qq.com/s/GpXovlWS_MAKdoRv3OAjCw
GDCdownload(query)
GDCprepare(query, save = T,save.filename = "TCGA-BLCA_SNP.Rdata")
# 
# library(TCGAbiolinks)
# # BiocManager::install("PoisonAlien/TCGAmutations")
# library(TCGAmutations)
# tmp=as.data.frame(tcga_available())
# 
# GBMmut <- GDCquery_Maf(tumor = "BLCA", pipelines = "mutect2")#下载突变数据

rm(list = ls())
load("TCGA-BLCA_SNP.Rdata")

file<-list.files(pattern = "Hs.gmt")
geneset<-data.frame()
for (i in 1:length(file)) {
  # i=1
  o<- clusterProfiler::read.gmt(file[i]) 
  geneset<-rbind(geneset,o)
}
gene<-unique(gene_set$gene)
gene_maf<-data[data$Hugo_Symbol%in%gene,]
length(unique(gene_maf$Hugo_Symbol))

library(maftools);maf = read.maf(maf = gene_maf)
vc_cols = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
pdf("mutation_genes.pdf",onefile = F,width = 10,height = 15)
oncoplot(maf = maf, colors = vc_cols, 
         top = 50,
         draw_titv = F)
dev.off()
