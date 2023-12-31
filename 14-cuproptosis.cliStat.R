
trainFile="risk.train.txt"    #train组数据文件
testFile="risk.test.txt"      #test组数据文件
cliFile="clinical.txt"        #临床数据文件
setwd("C:\\Users\\86188\\Desktop\\LUAD-Cu\\14.cliStat")      #设置工作目录

#读取train组数据文件
train=read.table(trainFile, header=T, sep="\t", check.names=F, row.names=1)
#读取test组数据文件
test=read.table(testFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床输入文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))

#合并数据
trainCli=cli[row.names(train),]
trainCli=cbind(trainCli, Type="Train")
testCli=cli[row.names(test),]
testCli=cbind(testCli, Type="Test")
rt=rbind(trainCli, testCli)

#对临床性状进行循环,得到train组和test组临床性状差异的pvalue
cliStatOut=data.frame()
for(i in 1:(ncol(rt)-1)){
  nameStat=colnames(rt)[i]
  tableStat=table(rt[,c(nameStat,"Type")])
  tableStatSum=cbind(Total=rowSums(tableStat), tableStat)
  tableStatRatio=prop.table(tableStatSum,2)
  tableStatRatio=round(tableStatRatio*100,2)
  tableStatPaste=paste(tableStatSum,"(",tableStatRatio,"%)",sep="")
  tableStatOut=matrix(tableStatPaste,ncol=3,dimnames=dimnames(tableStatSum))
  pStat=chisq.test(tableStat[row.names(tableStat)!="unknow",])
  pValueStat=round(pStat$p.value,4)
  pValueCol=c(pValueStat,rep(" ",(nrow(tableStatOut)-1)) )
  tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatOut),tableStatOut,Pvalue=pValueCol)
  cliStatOut=rbind(cliStatOut,tableStatOut)
}

#输出临床统计的结果
write.table(cliStatOut,file="cliStat.result.xls",sep="\t",quote=F,row.names=F)

