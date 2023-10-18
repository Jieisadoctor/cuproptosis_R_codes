#install.packages("pheatmap")

library(pheatmap)
setwd("C:\\Users\\86188\\Desktop\\LUAD-Cu\\15.riskPlot")             #设置工作目录

bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null,heatmapFile=null){
  rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)   #读取输入文件
  rt=rt[order(rt$riskScore),]                                            #按照riskScore对样品排序
  
  #绘制风险曲线
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=riskScoreFile,width = 10,height = 3.5)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("#336699",lowLength),rep("#996699",highLength)) )
  abline(h=median(rt$riskScore),v=lowLength,lty=2)
  legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("#996699","#336699"),cex=1.2)
  dev.off()
  
  #绘制生存状态图
  color=as.vector(rt$fustat)
  color[color==1]="#996699"
  color[color==0]="#336699"
  pdf(file=survStatFile,width = 10,height = 3.5)
  plot(rt$futime, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("#996699","#336699"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  
  #绘制风险热图
  rt1=rt[c(3:(ncol(rt)-2))]
  # rt1=log2(rt1+1)
  rt1=t(rt1)
  annotation=data.frame(type=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  pdf(file=heatmapFile,width = 10,height = 3)
  pheatmap(rt1, 
           annotation=annotation, 
           cluster_cols = FALSE,
           fontsize_row=11,
           show_colnames = F,
           fontsize_col=3,
           color = colorRampPalette(c("#336699", "white", "#996699"))(50) )
  dev.off()
}
bioRiskPlot(inputFile="risk.test.txt",riskScoreFile="geo.riskScore.pdf",survStatFile="geo.survStat.pdf",heatmapFile="geo.heatmap.pdf")
bioRiskPlot(inputFile="risk.train.txt",riskScoreFile="tcga.riskScore.pdf",survStatFile="tcga.survStat.pdf",heatmapFile="tcga.heatmap.pdf")

