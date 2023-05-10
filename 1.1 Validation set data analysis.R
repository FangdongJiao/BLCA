rm(list = ls())
myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");myfilepath
setwd(myfilepath)
#old-----------------------------------------------------------------
#GEO-GSE181507
{
  rm(list = ls())
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
  library(GEOquery)
  gse<-getGEO("GSE181507",getGPL = T)
  show(gse)
  save( gse, file = 'GSE181507_gset.Rdata' )
  
  load("GSE181507_gset.Rdata")
  a<-gse[[1]]
  show(a)
  expma<-exprs(a)
  id<-a@featureData@data
  symbol<-cbind(id$ID,id$SPOT_ID.1)
  symbol<-data.frame(symbol)
  n<-length(rownames(id))
  names(symbol)
  x<-strsplit(symbol$X2," // ")
  for (t in 1:n) {symbol$x4[t]<-x[[t]][3]}
  symbol1<-symbol[,-2]
  for (t in 1:nrow(symbol1)) {symbol1$x4[t]<-strsplit(symbol1$x4[t],"\\(")[[1]][2]}
  for (t in 1:nrow(symbol1)) {symbol1$x4[t]<-strsplit(symbol1$x4[t],"\\)")[[1]][1]}
  
  expma<-data.frame(expma)
  expma$X1<-rownames(expma)
  expr<-merge(symbol1,expma,by="X1")
  expr<-data.frame(expr)
  expr<-expr[,-1]
  mydf<-aggregate(expr[2:ncol(expr)],
                  list(expr$x4),
                  FUN = mean)
  rownames(mydf)<-mydf$Group.1
  GSE181507_matrix<-mydf[,-1]
  boxplot(GSE181507_matrix[1:50,])
  
  pdata<-pData(a)
  GSE181507_pdata<-pdata
  # group<-cbind(pdata$geo_accession,
  #              pdata$`age:ch1`,
  #              pdata$`disease_free_survival_event:ch1`,
  #              pdata$`tumor stage:ch1`)
  # group<-data.frame(group)
  # names(group)<-c("GSM","status","time","stage")
  # GSE181507_group<-group
  
  load("protein_coding.gtf.Rda")
  GSE181507_matrix<-GSE181507_matrix[protein_coding.gtf$gene_name,]
  GSE181507_matrix<-na.omit(GSE181507_matrix)
  
  save(GSE181507_matrix,GSE181507_pdata,file = "GSE181507.Rdata")
  rm(list = ls())
  load("GSE181507.Rdata")
}

#GEO-GSE13507
{
  rm(list = ls())
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
  library(GEOquery)
  gse<-getGEO("GSE13507",getGPL = T)
  show(gse)
  save( gse, file = 'GSE13507_gset.Rdata' )
  
  load("GSE13507_gset.Rdata")
  a<-gse[[1]]
  show(a)
  expma<-exprs(a)
  id<-gse$GSE13507_series_matrix.txt.gz@featureData@data
  symbol<-cbind(id$ID,id$Symbol)
  symbol<-data.frame(symbol)
  # n<-length(rownames(id))
  # x<-strsplit(symbol$X2," // ")
  # for (t in 1:n) {symbol$x4[t]<-x[[t]][2]}
  symbol1<-symbol
  
  expma<-data.frame(expma)
  expma$X1<-rownames(expma)
  expr<-merge(symbol1,expma,by="X1")
  expr<-data.frame(expr)
  expr<-expr[,-1]
  mydf<-aggregate(expr[2:ncol(expr)],
                  list(expr$X2),
                  FUN = mean)
  rownames(mydf)<-mydf$Group.1
  GSE13507_matrix<-mydf[,-1]
  boxplot(GSE13507_matrix[1:50,1:50])
  
  pdata<-pData(a)
  GSE13507_pdata<-pdata
  # group<-cbind(pdata$geo_accession,
  #              pdata$`60 months overall survival, status:ch1`,
  #              pdata$`60 months overall survival, time:ch1`,
  #              pdata$`gender:ch1`,
  #              pdata$`tumor number:ch1`)
  
  load("protein_coding.gtf.Rda")
  GSE13507_matrix<-GSE13507_matrix[protein_coding.gtf$gene_name,]
  GSE13507_matrix<-na.omit(GSE13507_matrix)
  
  save(GSE13507_matrix,GSE13507_pdata,file = "GSE13507.Rdata")
  rm(list = ls())
  load("GSE13507.Rdata")
}

#GEO-GSE48276
{
  rm(list = ls())
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
  library(GEOquery)
  gse<-getGEO("GSE48276",getGPL = T)
  show(gse)
  save( gse, file = 'GSE48276_gset.Rdata' )
  
  load("GSE48276_gset.Rdata")
  a<-gse[[1]]
  show(a)
  expma<-exprs(a)
  id<-gse$GSE48276_series_matrix.txt.gz@featureData@data
  symbol<-cbind(id$ID,id$Symbol)
  symbol<-data.frame(symbol)
  symbol1<-symbol
  
  expma<-data.frame(expma)
  expma$X1<-rownames(expma)
  expr<-merge(symbol1,expma,by="X1")
  expr<-data.frame(expr)
  expr<-expr[,-1]
  mydf<-aggregate(expr[2:ncol(expr)],
                  list(expr$X2),
                  FUN = mean)
  rownames(mydf)<-mydf$Group.1
  GSE48276_matrix<-mydf[,-1]
  boxplot(GSE48276_matrix[1:50,1:50])
  
  pdata<-pData(a)
  GSE48276_pdata<-pdata
  
  load("protein_coding.gtf.Rda")
  GSE48276_matrix<-GSE48276_matrix[protein_coding.gtf$gene_name,]
  GSE48276_matrix<-na.omit(GSE48276_matrix)
  
  save(GSE48276_matrix,GSE48276_pdata,file = "GSE48276.Rdata")
  rm(list = ls())
  load("GSE48276.Rdata")
}

#GEO-GSE19915
{
  rm(list = ls())
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
  library(GEOquery)
  gse<-getGEO("GSE19915",getGPL = T)
  show(gse)
  save( gse, file = 'GSE19915_gset.Rdata' )
  
  a<-gse[[1]]
  show(a)
  expma<-exprs(a)
  id<-a@featureData@data
  symbol<-cbind(id$ID,id$GeneSymbol)
  symbol<-data.frame(symbol)
  n<-length(rownames(id))
  symbol$X2
  x<-strsplit(symbol$X2,"\\|")
  for (t in 1:n) {symbol$x4[t]<-x[[t]][1]}
  symbol1<-symbol[,c(1,3)]
  
  expma<-data.frame(expma)
  expma$X1<-rownames(expma)
  expr<-merge(symbol1,expma,by="X1")
  expr<-data.frame(expr)
  expr<-expr[,-1]
  mydf<-aggregate(expr[2:ncol(expr)],
                  list(expr$x4),
                  FUN = mean)
  rownames(mydf)<-mydf$Group.1
  GSE19915_matrix<-mydf[,-1]
  boxplot(GSE19915_matrix[1:50,])
  
  # load("protein_coding.gtf.Rda")
  # GSE19915_matrix<-GSE19915_matrix[protein_coding.gtf$gene_name,]
  # GSE19915_matrix<-na.omit(GSE19915_matrix)
  
  pdata<-pData(a)
  GSE19915_pdata<-pdata
  
  
  save(GSE19915_matrix,GSE19915_pdata,file = "GSE19915.Rdata")
  rm(list = ls())
  load("GSE19915.Rdata")
}

#new-----------------------------------------------------------------
#GEO-GSE37815
{
  rm(list = ls())
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
  library(GEOquery)
  gse<-getGEO("GSE37815",getGPL = F)
  show(gse)
  save( gse, file = 'GSE37815_gset.Rdata' )
  
  load("GSE37815_gset.Rdata")
  # GSE37815<-gse
  # 
  # load("GSE13507_gset.Rdata")
  # id<-gse$GSE13507_series_matrix.txt.gz@featureData@data
  
  rm(list = ls())
  library(data.table)
  exp<-data.frame(fread("GSE37815_non-normalized.txt.gz"))
  exp<-exp[-1,]
  id<-fread("GPL6102-11574.txt")
  symbol<-cbind(id$ID,id$Symbol)
  symbol<-data.frame(symbol)
  # n<-length(rownames(id))
  # x<-strsplit(symbol$X2," // ")
  # for (t in 1:n) {symbol$x4[t]<-x[[t]][2]}
  symbol1<-symbol
  
  expma<-data.frame(exp)
  names(expma)[1]<-"X1"
  expr<-merge(symbol1,expma,by="X1")
  expr<-data.frame(expr)
  expr<-expr[,-1]
  mydf<-expr[!duplicated(expr$X2),]
  rownames(mydf)<-mydf[,1]
  GSE37815_matrix<-mydf[,-1]
  
  GSE37815_matrix<-apply(GSE37815_matrix, 2, as.numeric)
  rownames(GSE37815_matrix)<-mydf$X2
  
  boxplot(GSE37815_matrix[1:50,])
  d<-log2(GSE37815_matrix+0.00001)
  boxplot(d[1:50,])
  GSE37815_matrix<-d
  
  a<-gse[[1]]
  pdata<-pData(a)
  GSE37815_pdata<-pdata
  
  
  save(GSE37815_matrix,GSE37815_pdata,file = "GSE37815.Rdata")
  rm(list = ls())
  load("GSE37815.Rdata")
}

#GEO-GSE31684×
{
  rm(list = ls())
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
  library(GEOquery)
  gse<-getGEO("GSE31684",getGPL = T)
  show(gse)
  save( gse, file = 'GSE31684_gset.Rdata' )
  load("GSE31684_gset.Rdata")
  a<-gse[[1]]
  show(a)
  expma<-exprs(a)
  library(data.table)
  file<-list.files("./GSE31684_RAW")
  for (i in 1:length(file)) {
    i=1
    expression<-fread(paste("./GSE31684_RAW/",file[i],sep = ""))

  }

  expma
  
  
  
  id<-fread("GPL570-55999.txt")
  symbol<-cbind(id$ID,id$`Gene Symbol`)
  symbol<-data.frame(symbol)
  n<-length(rownames(id))
  head(symbol$X2)
  x<-strsplit(symbol$X2,"///")
  for (t in 1:n) {symbol$x4[t]<-x[[t]][1]}
  symbol1<-symbol[,c(1,3)]
  
  expma<-data.frame(expma)
  expma$X1<-rownames(expma)
  expr<-merge(symbol1,expma,by="X1")
  expr<-data.frame(expr)
  expr<-expr[,-1]
  mydf<-aggregate(expr[2:ncol(expr)],
                  list(expr$x4),
                  FUN = mean)
  rownames(mydf)<-mydf$Group.1
  GSE31684_matrix<-mydf[,-1]
  boxplot(GSE31684_matrix[1:50,])
  
  # load("protein_coding.gtf.Rda")
  # GSE31684_matrix<-GSE31684_matrix[protein_coding.gtf$gene_name,]
  # GSE31684_matrix<-na.omit(GSE31684_matrix)
  
  pdata<-pData(a)
  GSE31684_pdata<-pdata
  
  
  save(GSE31684_matrix,GSE31684_pdata,file = "GSE31684.Rdata")
  rm(list = ls())
  load("GSE31684.Rdata")
}

#GEO-GSE19423
{
  rm(list = ls())
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
  library(GEOquery)
  gse<-getGEO("GSE19423",getGPL = T)
  show(gse)
  save( gse, file = 'GSE19423_gset.Rdata' )
  
  load("GSE19423_gset.Rdata")
  library(data.table)
  exp<-data.frame(fread("GSE19423_raw_data.txt.gz"))
  id<-fread("GPL6102-11574.txt")
  symbol<-cbind(id$ID,id$Symbol)
  symbol<-data.frame(symbol)
  # n<-length(rownames(id))
  # x<-strsplit(symbol$X2," // ")
  # for (t in 1:n) {symbol$x4[t]<-x[[t]][2]}
  symbol1<-symbol
  
  expma<-data.frame(exp)
  names(expma)[1]<-"X1"
  expr<-merge(symbol1,expma,by="X1")
  expr<-data.frame(expr)
  expr<-expr[,-1]
  mydf<-expr[!duplicated(expr$X2),]
  rownames(mydf)<-mydf$Group.1
  GSE19423_matrix<-mydf[,-1]
  boxplot(GSE19423_matrix[1:50,])
  
  boxplot(GSE19423_matrix[1:50,])
  d<-log2(GSE19423_matrix+0.00001)
  boxplot(d[1:50,])
  GSE19423_matrix<-d
  
  a<-gse[[1]]
  pdata<-pData(a)
  GSE19423_pdata<-pdata
  
  save(GSE19423_matrix,GSE19423_pdata,file = "GSE19423.Rdata")
  rm(list = ls())
  load("GSE19423.Rdata")
}


#GEO-GSE149582
{
  rm(list = ls())
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
  library(GEOquery)
  gse<-getGEO("GSE149582",getGPL = T)
  show(gse)
  save( gse, file = 'GSE149582_gset.Rdata' )
  
  load("GSE149582_gset.Rdata")
  a<-gse[[1]]
  show(a)
  expma<-exprs(a)
  id<-gse$GSE149582_series_matrix.txt.gz@featureData@data
  symbol<-cbind(id$ID,id$Symbol)
  symbol<-data.frame(symbol)
  # n<-length(rownames(id))
  # x<-strsplit(symbol$X2," // ")
  # for (t in 1:n) {symbol$x4[t]<-x[[t]][2]}
  symbol1<-symbol
  
  expma<-data.frame(expma)
  expma$X1<-rownames(expma)
  expr<-merge(symbol1,expma,by="X1")
  expr<-data.frame(expr)
  expr<-expr[,-1]
  mydf<-aggregate(expr[2:ncol(expr)],
                  list(expr$X2),
                  FUN = mean)
  rownames(mydf)<-mydf$Group.1
  GSE149582_matrix<-mydf[,-1]
  boxplot(GSE149582_matrix[1:50,1:4])
  
  pdata<-pData(a)
  GSE149582_pdata<-pdata
  # group<-cbind(pdata$geo_accession,
  #              pdata$`60 months overall survival, status:ch1`,
  #              pdata$`60 months overall survival, time:ch1`,
  #              pdata$`gender:ch1`,
  #              pdata$`tumor number:ch1`)
  
  save(GSE149582_matrix,GSE149582_pdata,file = "GSE149582.Rdata")
  rm(list = ls())
  load("GSE149582.Rdata")
}


#GEO-GSE39280
{
  rm(list = ls())
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
  library(GEOquery)
  gse<-getGEO("GSE39280",getGPL = T)
  show(gse)
  save( gse, file = 'GSE39280_gset.Rdata' )
  
  load("GSE39280_gset.Rdata")
  a<-gse[[1]]
  show(a)
  expma<-exprs(a)
  id<-gse$GSE39280_series_matrix.txt.gz@featureData@data
  symbol<-cbind(id$ID,id$Gene)
  symbol<-data.frame(symbol)
  # n<-length(rownames(id))
  # x<-strsplit(symbol$X2," // ")
  # for (t in 1:n) {symbol$x4[t]<-x[[t]][2]}
  symbol1<-symbol
  
  expma<-data.frame(expma)
  expma$X1<-rownames(expma)
  expr<-merge(symbol1,expma,by="X1")
  expr<-data.frame(expr)
  expr<-expr[,-1]
  mydf<-aggregate(expr[2:ncol(expr)],
                  list(expr$X2),
                  FUN = mean)
  rownames(mydf)<-mydf$Group.1
  GSE39280_matrix<-mydf[,-1]
  boxplot(GSE39280_matrix[1:50,1:4])
  
  pdata<-pData(a)
  GSE39280_pdata<-pdata
  # group<-cbind(pdata$geo_accession,
  #              pdata$`60 months overall survival, status:ch1`,
  #              pdata$`60 months overall survival, time:ch1`,
  #              pdata$`gender:ch1`,
  #              pdata$`tumor number:ch1`)
  
  save(GSE39280_matrix,GSE39280_pdata,file = "GSE39280.Rdata")
  rm(list = ls())
  load("GSE39280.Rdata")
}

#GEO-GSE39281 ✔
{
  rm(list = ls())
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
  library(GEOquery)
  gse<-getGEO("GSE39281",getGPL = T)
  show(gse)
  save( gse, file = 'GSE39281_gset.Rdata' )
  
  load("GSE39281_gset.Rdata")
  a<-gse[[1]]
  show(a)
  expma<-exprs(a)
  id<-a@featureData@data
  symbol<-cbind(id$ID,id$GENE_SYMBOL)
  symbol<-data.frame(symbol)
  n<-length(rownames(id))
  symbol$X2
  x<-strsplit(symbol$X2,"\\|")
  for (t in 1:n) {symbol$x4[t]<-x[[t]][1]}
  symbol1<-symbol[,c(1,3)]
  
  expma<-data.frame(expma)
  expma$X1<-rownames(expma)
  expr<-merge(symbol1,expma,by="X1")
  expr<-data.frame(expr)
  expr<-expr[,-1]
  mydf<-aggregate(expr[2:ncol(expr)],
                  list(expr$x4),
                  FUN = mean)
  rownames(mydf)<-mydf$Group.1
  GSE39281_matrix<-mydf[,-1]
  boxplot(GSE39281_matrix[1:50,])
  
  # load("protein_coding.gtf.Rda")
  # GSE39281_matrix<-GSE39281_matrix[protein_coding.gtf$gene_name,]
  # GSE39281_matrix<-na.omit(GSE39281_matrix)
  
  pdata<-pData(a)
  GSE39281_pdata<-pdata
  
  
  save(GSE39281_matrix,GSE39281_pdata,file = "GSE39281.Rdata")
  rm(list = ls())
  load("GSE39281.Rdata")
}

#GEO-GSE48075
{
  rm(list = ls())
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
  library(GEOquery)
  gse<-getGEO("GSE48075",getGPL = T)
  show(gse)
  save( gse, file = 'GSE48075_gset.Rdata' )
  
  load("GSE48075_gset.Rdata")
  a<-gse[[1]]
  show(a)
  expma<-exprs(a)
  id<-a@featureData@data
  symbol<-cbind(id$ID,id$Symbol)
  symbol<-data.frame(symbol)
  n<-length(rownames(id))
  head(symbol)
  # x<-strsplit(symbol$X2," // ")
  # for (t in 1:n) {symbol$x4[t]<-x[[t]][3]}
  # symbol1<-symbol[,-2]
  symbol1<-symbol
  expma<-data.frame(expma)
  expma$X1<-rownames(expma)
  expr<-merge(symbol1,expma,by="X1")
  expr<-data.frame(expr)
  expr<-expr[,-1]
  mydf<-expr[!duplicated(expr$X2),]
  rownames(mydf)<-mydf[,1]
  GSE48075_matrix<-mydf[,-1]
  boxplot(GSE48075_matrix[1:50,])
  
  pdata<-pData(a)
  GSE48075_pdata<-pdata
  # group<-cbind(pdata$geo_accession,
  #              pdata$`age:ch1`,
  #              pdata$`disease_free_survival_event:ch1`,
  #              pdata$`tumor stage:ch1`)
  # group<-data.frame(group)
  # names(group)<-c("GSM","status","time","stage")
  # GSE48075_group<-group
  save(GSE48075_matrix,GSE48075_pdata,file = "GSE48075.Rdata")
  rm(list = ls())
  load("GSE48075.Rdata")
}

#GEO-GSE69795
{
  rm(list = ls())
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
  library(GEOquery)
  gse<-getGEO("GSE69795",getGPL = T)
  show(gse)
  save( gse, file = 'GSE69795_gset.Rdata' )
  
  load("GSE69795_gset.Rdata")
  a<-gse[[1]]
  show(a)
  expma<-exprs(a)
  id<-gse$GSE69795_series_matrix.txt.gz@featureData@data
  symbol<-cbind(id$ID,id$Symbol)
  symbol<-data.frame(symbol)
  # n<-length(rownames(id))
  # x<-strsplit(symbol$X2," // ")
  # for (t in 1:n) {symbol$x4[t]<-x[[t]][2]}
  symbol1<-symbol
  
  expma<-data.frame(expma)
  expma$X1<-rownames(expma)
  expr<-merge(symbol1,expma,by="X1")
  expr<-data.frame(expr)
  expr<-expr[,-1]
  mydf<-aggregate(expr[2:ncol(expr)],
                  list(expr$X2),
                  FUN = mean)
  rownames(mydf)<-mydf$Group.1
  GSE69795_matrix<-mydf[,-1]
  boxplot(GSE69795_matrix[1:50,1:50])
  
  pdata<-pData(a)
  GSE69795_pdata<-pdata
  # group<-cbind(pdata$geo_accession,
  #              pdata$`60 months overall survival, status:ch1`,
  #              pdata$`60 months overall survival, time:ch1`,
  #              pdata$`gender:ch1`,
  #              pdata$`tumor number:ch1`)
  
  save(GSE69795_matrix,GSE69795_pdata,file = "GSE69795.Rdata")
  rm(list = ls())
  load("GSE69795.Rdata")
}


#GEO-GSE37816
{
  rm(list = ls())
  Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
  library(GEOquery)
  gse<-getGEO("GSE37816",getGPL = T)
  show(gse)
  save( gse, file = 'GSE37816_gset.Rdata' )
  
  load("GSE37816_gset.Rdata")
  a<-gse[[1]]
  show(a)
  expma<-exprs(a)
  id<-gse$GSE37816_series_matrix.txt.gz@featureData@data
  symbol<-cbind(id$ID,id$Symbol)
  symbol<-data.frame(symbol)
  # n<-length(rownames(id))
  # x<-strsplit(symbol$X2," // ")
  # for (t in 1:n) {symbol$x4[t]<-x[[t]][2]}
  symbol1<-symbol
  
  expma<-data.frame(expma)
  expma$X1<-rownames(expma)
  expr<-merge(symbol1,expma,by="X1")
  expr<-data.frame(expr)
  expr<-expr[,-1]
  mydf<-aggregate(expr[2:ncol(expr)],
                  list(expr$X2),
                  FUN = mean)
  rownames(mydf)<-mydf$Group.1
  GSE37816_matrix<-mydf[,-1]
  boxplot(GSE37816_matrix[1:50,1:4])
  
  pdata<-pData(a)
  GSE37816_pdata<-pdata
  # group<-cbind(pdata$geo_accession,
  #              pdata$`60 months overall survival, status:ch1`,
  #              pdata$`60 months overall survival, time:ch1`,
  #              pdata$`gender:ch1`,
  #              pdata$`tumor number:ch1`)
  
  save(GSE37816_matrix,GSE37816_pdata,file = "GSE37816.Rdata")
  rm(list = ls())
  load("GSE37816.Rdata")
}