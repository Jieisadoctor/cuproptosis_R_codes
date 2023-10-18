
library(limma)
setwd("C:\\Users\\86188\\Desktop\\LUAD-Cu\\01.download")          #设置工作目录

#读取symbol文件,并对数据进行处理
rt = read.table("symbol.txt",header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
lncRNA=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
lncRNA=avereps(lncRNA)
lncRNA=lncRNA[rowMeans(lncRNA)>0.5,]
lncRNA=log2(lncRNA+1)
lncRNA=normalizeBetweenArrays(lncRNA)
lncRNAOUT=cbind(id=rownames(lncRNA),lncRNA)
write.table(lncRNAOUT,file="symbol_correct.txt",sep="\t",quote=F,row.names=F)

