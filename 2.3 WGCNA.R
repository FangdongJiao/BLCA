myfilepath<-dirname(rstudioapi::getActiveDocumentContext()$path);myfilepath
str<-strsplit(myfilepath,"/R-code");myfilepath<-paste(str,"/LZX文件",sep = "");
setwd(myfilepath)

rm(list = ls())
load("TCGA-BLCA.Rdata")
table(TCGA_group$group_list)
tumor<-TCGA_group[TCGA_group$group_list=="tumor",]

DEG_exp<-TCGA_exp[,tumor$GSM]
library(WGCNA)
options(stringsAsFactors = FALSE)
# 指允许R语言程序最大线程运行
allowWGCNAThreads()
samples=TCGA_group[,c(2,1)]
datExpr=as.data.frame(t(DEG_exp))
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##软阈值筛选##
powers = c(seq(1,10,by = 1), seq(12, 20, by = 2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("β.pdf",width=8,height=6,onefile = FALSE) 
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="green");
abline(h=0.8,col="green")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="green")
dev.off()
sft$powerEstimate
dev.off()

k <- softConnectivity(datE=datExpr,power=5) 
pdf("K.pdf",width=8,height=6,onefile = FALSE) 
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()

net = blockwiseModules(datExpr, power = 5, maxBlockSize = 6000,
                       TOMType = "unsigned", minModuleSize = 500,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       verbose = 3)
table(net$colors)
# cor<-stats::cor
#sft$powerEstimate
mergedColors = labels2colors(net$colors)

pdf("modele.pdf",width=10,height=6,onefile = FALSE) 
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()

moduleLabels = net$colors;moduleLabels
moduleColors = labels2colors(net$colors);moduleColors
table(moduleColors)
MEs_col = net$MEs

load("ssGSEA.Rdata")
k=data.frame(t(gsva_matrix))

moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEs = orderMEs(MEs0)
MEsWW<-MEs

k2<-k[rownames(MEsWW),]
modTraitCor = cor(MEsWW,k2,use = "p")

colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)

textMatrix = paste(signif(modTraitCor, 2), 
                   "\n(", signif(modTraitP, 1), ")", 
                   sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("WGCNA_COR.pdf",width=10,height=10,onefile = FALSE) 
par(mar=c(20,20,20,20))
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(k2), 
               yLabels = names(MEsWW), cex.lab = 0.9,  yColorWidth=0.01, 
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix,
               setStdMargins = FALSE, cex.text = 0.5, 
               zlim = c(-1,1), 
               main = paste("Module-trait relationships"))

dev.off()


f1<-c("yellow","turquoise","green","red","blue")#1
f1
num1<-c()
GENE<-list()
for (i in 1:length(f1)) {
  # i=1
  names(data.frame(modTraitCor))
  module = f1[i]
  moduleGenes = moduleColors == module
  Meta_gene=data.frame(moduleGenes)
  print(table(Meta_gene))
  g=data.frame(colnames(datExpr))
  gene<-g[Meta_gene==T,]
  gene<-data.frame(gene)
  
  GENE[[i]]<-gene
  names(GENE)[i]<-f1[i]
  
  num1[i]<-length(gene$gene)
  names(num1)[i]<-f1[i]
  write.table(gene,paste("Pathway_",f1[i],".xls",sep = ""),quote = F,sep = "\t",row.names = F)
}
save(GENE,file = "gene_modelu.Rda")