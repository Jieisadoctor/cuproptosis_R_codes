#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
# install.packages("ggExtra")


#引用包
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

expFile="RiskScore.txt"     #表达数据文件
tmbFile="TMB.txt"         #肿瘤突变符合文件
setwd("C:\\Users\\86188\\Desktop\\LUAD-Cu\\27.TMBcor")     #设置工作目录
#读取肿瘤突变负荷文件
tmb=read.table("TMB.txt", header=T, sep="\t", check.names=F, row.names=1)

#读取风险数据文件
risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(tmb, risk)

#相关性分析
x=as.numeric(data[,"riskScore"])
y=log2(as.numeric(data[,"TMB"])+1)
df1=as.data.frame(cbind(x,y))
corT=cor.test(x, y, method="spearman")
p1=ggplot(df1, aes(x, y)) + 
  xlab(paste0("Risk score"))+ylab("Tumor mutation burden")+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))

#输出相关性图形
pdf(file="cor.pdf",width=4,height=3.8)
print(p1)
dev.off()
