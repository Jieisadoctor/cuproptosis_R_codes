#install.packages("survival")
#install.packages("regplot")
#install.packages("rms")


#引用包
library(survival)
library(regplot)
library(rms)

riskFile="risk.all.txt"      #风险文件
riskFiletrain="risk.train.txt"      #风险文件
riskFiletest="risk.test.txt"      #风险文件
cliFile="clinical.txt"       #临床数据文件
setwd("C:\\Users\\86188\\Desktop\\LUAD-Cu\\17.Nomo")     #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1[,c("futime", "fustat", "riskScore")], cli)

#训练集数据
risktrain=read.table(riskFiletrain, header=T, sep="\t", check.names=F, row.names=1)
clitrain=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
clitrain=clitrain[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
clitrain$Age=as.numeric(clitrain$Age)
samSampletrain=intersect(row.names(risktrain), row.names(clitrain))
risk1train=risktrain[samSampletrain,,drop=F]
clitrain=clitrain[samSampletrain,,drop=F]
rttrain=cbind(risk1train[,c("futime", "fustat", "riskScore")], clitrain)

#测试集数据
risktest=read.table(riskFiletest, header=T, sep="\t", check.names=F, row.names=1)
clitest=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
clitest=clitest[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
clitest$Age=as.numeric(clitest$Age)
samSampletest=intersect(row.names(risktest), row.names(clitest))
risk1test=risktest[samSampletest,,drop=F]
clitest=clitest[samSampletest,,drop=F]
rttest=cbind(risk1test[,c("futime", "fustat", "riskScore")], clitest)

#绘制列线图
res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
             plots = c("density", "boxes"),
             clickable=F,
             title="",
             points=TRUE,
             droplines=TRUE,
             observation=rt[20,],
             rank="sd",
             failtime = c(1,3,5),
             prfail = F)

#列线图风险打分
nomoRisk=predict(res.cox, newdata=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)

#校准曲线
pdf(file="calibration.pdf", width=5, height=5)
#1年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/2), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)
#3年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/2), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)
#5年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/2), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()

#训练集数据
#列线图风险打分
nomoRisktrain=predict(res.cox, newdata=rttrain, type="risk")
rttrain=cbind(risk1train, Nomogram=nomoRisktrain)
outTabtrain=rbind(ID=colnames(rttrain), rttrain)
write.table(outTabtrain, file="nomoRisktrain.txt", sep="\t", col.names=F, quote=F)

#校准曲线
pdf(file="calibrationtrain.pdf", width=5, height=5)
#1年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rttrain, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rttrain)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)
#3年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rttrain, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rttrain)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)
#5年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rttrain, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rttrain)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()

#测试集数据
#列线图风险打分
nomoRisktest=predict(res.cox, newdata=rttest, type="risk")
rttest=cbind(risk1test, Nomogram=nomoRisktest)
outTabtest=rbind(ID=colnames(rttest), rttest)
write.table(outTabtest, file="nomoRisktest.txt", sep="\t", col.names=F, quote=F)

#校准曲线
pdf(file="calibrationtest.pdf", width=5, height=5)
#1年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rttest, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rttest)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)
#3年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rttest, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rttest)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)
#5年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rttest, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rttest)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()
