
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")

logFoldChange=0.5               #logFC过滤阈值
adjustP=0.05                   #矫正后p值阈值
conNum=54                     #control组样品数目
treatNum=497                  #treat组样品数目

library(limma)
library(pheatmap)
setwd("C:\\Users\\86188\\Desktop\\LUAD-Cu\\10.lncRNAdiff")        #设置工作目录
rt=read.table("lncRNA.txt",sep="\t",header=T,check.names=F)    #读取输入文件
#数据处理，如果一个基因有多行，取均值

rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
#数据矫正
# rt=normalizeBetweenArrays(as.matrix(rt))
#如果下载的数据没有取log,需要把下面这行命令前面的#号去掉
# rt=log2(rt+1)

#读取铜死亡lncRNA
lncRNA=read.table("immuneLncRNAexp.txt",sep="\t",header=T,check.names=F)
lncRNA_read=lncRNA[,1]
rt=rt[lncRNA_read,]

#differential
modType=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(modType))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="mirnaAll.xls",sep="\t",quote=F)

#write table
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="mirnaDiff.xls",sep="\t",quote=F,col.names=F)
write.table(diffSigOut,file="mirnaDiff.txt",sep="\t",quote=F,col.names=F)
DiffexpOUT=rt[row.names(diffSig),]
DiffexpOUT=cbind(id=rownames(DiffexpOUT),DiffexpOUT)
write.table(DiffexpOUT,file="DiffexpOUT.txt",sep="\t",quote=F,row.names=F)

#输出差异lncRNA
difflncRNA=as.matrix(lncRNA)
rownames(lncRNA)=lncRNA[,1]
exp2=lncRNA[,2:ncol(lncRNA)]
dimnames2=list(rownames(exp2),colnames(exp2))
difflncRNA=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
difflncRNAexp <- difflncRNA[rownames(diffSig),]
difflncRNAexp_out=rbind(id=colnames(difflncRNAexp),difflncRNAexp)
write.table(difflncRNAexp_out,file="difflncRNAexp.txt",sep="\t",quote=F,col.names=F)

#绘制火山图
pdf(file="mirnaVol.pdf",height=5,width=5)
xMax=max(abs(allDiff$logFC))
yMax=max(-log10(allDiff$adj.P.Val))
plot(allDiff$logFC, -log10(allDiff$adj.P.Val), xlab="logFC",ylab="-log10(adj.P.Val)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(allDiff, adj.P.Val<adjustP & logFC>logFoldChange)
points(diffSub$logFC, -log10(diffSub$adj.P.Val), pch=20, col="red",cex=0.8)
diffSub=subset(allDiff, adj.P.Val<adjustP & logFC<(-logFoldChange))
points(diffSub$logFC, -log10(diffSub$adj.P.Val), pch=20, col="blue",cex=0.8)
abline(v=0,lty=2,lwd=3)
dev.off()

#绘制差异基因热图ry(pheatmap)
hmExp=rt[rownames(diffSig),]
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(rt)
Type=as.data.frame(Type)
pdf(file="mirnaHeatmap.pdf",height=4,width=6)
pheatmap(hmExp, 
         annotation=Type, 
         # annotation_colors = list(Type = c(T = "#996699", N = "darkblue")),
         cluster_col = FALSE,
         color = colorRampPalette(c(rep("#336699",2), "white", rep("#996699",2)))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 12,
         fontsize_row=4,
         fontsize_col=10)
dev.off()

