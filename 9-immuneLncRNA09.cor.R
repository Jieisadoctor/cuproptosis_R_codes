#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library(limma)
setwd("C:\\Users\\86188\\Desktop\\LUAD-Cu\\09.immuneLncRNA")          #设置工作目录

corFilter=0.3             #相关系数过滤标准
pvalueFilter=0.05         #p值过滤标准

#读取lncRNA表达文件,并对数据进行处理
rt = read.table("lncRNA.txt",header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
lncRNA=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
lncRNA=avereps(lncRNA)
lncRNA=lncRNA[rowMeans(lncRNA)>0.5,]
group=sapply(strsplit(colnames(lncRNA),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
lncRNA=lncRNA[,group==0]

#读取免疫基因表达文件,并对数据进行处理
rt = read.table("ARGexp.txt",header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
immuneGene=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
immuneGene=avereps(immuneGene)
immuneGene=immuneGene[rowMeans(immuneGene)>0.5,]
group=sapply(strsplit(colnames(immuneGene),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
immuneGene=immuneGene[,group==0]

#相关性检验
outTab=data.frame()
for(i in row.names(lncRNA)){
  if(sd(lncRNA[i,])>0.5){
    for(j in row.names(immuneGene)){
      x=as.numeric(lncRNA[i,])
      y=as.numeric(immuneGene[j,])
      corT=cor.test(x,y)
      cor=corT$estimate
      pvalue=corT$p.value
      if((cor>corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(immuneGene=j,lncRNA=i,cor,pvalue,Regulation="postive"))
      }
      if((cor< -corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(immuneGene=j,lncRNA=i,cor,pvalue,Regulation="negative"))
      }
    }
  }
}
write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F)        #输出相关性结果
immuneLncRNA=unique(outTab[,"lncRNA"])
immuneLncRNAexp=lncRNA[immuneLncRNA,]
immuneLncRNAexp=rbind(ID=colnames(immuneLncRNAexp),immuneLncRNAexp)
write.table(immuneLncRNAexp,file="immuneLncRNAexp.txt",sep="\t",quote=F,col.names=F)

