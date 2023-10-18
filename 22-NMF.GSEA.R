
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

expFile="symbol.txt"             #表达数据文件
riskFile="risk.all.txt"      #风险文件
gmtFile="h.all.v7.5.1.symbols.gmt"     #基因集文件
setwd("C:\\Users\\86188\\Desktop\\LUAD-Cu\\22.GSEA")      #设置工作目录

#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))

#读取风险文件
Risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
data=data[,row.names(Risk)]

#高低风险比较，得到logFC
dataL=data[,row.names(Risk[Risk[,"risk"]=="low",])]
dataH=data[,row.names(Risk[Risk[,"risk"]=="high",])]
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
logFC=sort(logFC,decreasing=T)
genes=names(logFC)

#读入基因集文件
gmt=read.gmt(gmtFile)

#富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.25 & abs(kkTab$NES) > 1,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)

#输出高风险富集的图形
kkUp=kkTab[kkTab$NES>0,]
set.seed(122)
Ran1=sample(1:12, 6, replace=F)
kkUp1=kkUp[Ran1,]
kkUp2=kkUp[-Ran1,]
  showTerm=row.names(kkUp1)[1:nrow(kkUp1)]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in high risk group")
  pdf(file="GSEA.highRisk1.pdf", width=7, height=5.5)
  print(gseaplot)
  dev.off()
  
  showTerm=row.names(kkUp2)[1:nrow(kkUp2)]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in high risk group")
  pdf(file="GSEA.highRisk2.pdf", width=7, height=5.5)
  print(gseaplot)
  dev.off()
  

#输出低风险富集的图形
# termNum=5     #展示前5个通路
kkDown=kkTab[kkTab$NES<0,]
set.seed(123)
Ran1=sample(1:7, 4, replace=F)
kkDown1=kkDown[Ran1,]
kkDown2=kkDown[-Ran1,]
# if(nrow(kkDown)>=termNum){
  showTerm=row.names(kkDown1)[1:nrow(kkDown1)]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in low risk group")
  pdf(file="GSEA.lowRisk1.pdf", width=7, height=5.5)
  print(gseaplot)
  dev.off()
  
  showTerm=row.names(kkDown2)[1:nrow(kkDown2)]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in low risk group")
  pdf(file="GSEA.lowRisk2.pdf", width=7, height=5.5)
  print(gseaplot)
  dev.off()
# }


