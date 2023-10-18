
library(dplyr)
library(ggplot2)
setwd("C:\\Users\\86188\\Desktop\\LUAD-Cu\\17.3.clinical")    #设置工作目录
dat=read.table("clinical.txt",header=T,sep="\t",check.names=F,row.names = 1) 
risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)
dat[,"Age"]=ifelse(dat[,"Age"]=="unknow", "unknow", ifelse(dat[,"Age"]>65,">65","<=65"))

#样本合并
sameSample=intersect(row.names(dat),row.names(risk))
risk=risk[sameSample,]
dat=dat[sameSample,]
dat=cbind(dat,risk=risk[,"risk"])

# 按Risk分成High和Low，计算各列数值。
gname <- "risk"
vname <- setdiff(colnames(dat), gname)
pie.high <- pie.low <- list()
fisher.p <- c()

for (i in vname) {
  tmp <- table(dat[,gname], dat[,i])
  p <- format(fisher.test(tmp)$p.value,digits = 2)
  names(p) <- i
  fisher.p <- c(fisher.p, p)
  
  pie.dat <- 
    tmp %>% as.data.frame() %>% group_by(Var1) %>% mutate(Pct = Freq/sum(Freq)) %>% as.data.frame()
  
  # 表格内的两行对应Risk的两类：Risk high和Risk low
  pie.high[[i]] <- pie.dat[which(pie.dat$Var1 == "high"),]
  pie.low[[i]] <- pie.dat[which(pie.dat$Var1 == "low"),]
}
# 
# rthigh=rt[rt$risk=="high",]
# rthigh=rthigh[,1:ncol(rthigh)-1]
# rtlow=rt[rt$risk=="low",]
# rtlow=rtlow[,1:ncol(rtlow)-1]

#定义颜色
black  <- "#1E1E1B"
blue   <- "#3C4E98"
yellow <- "#E4DB36"
orange <- "#E19143"
green  <- "#57A12B"
cherry <- "#8D3A86"
red <- "#ff3333"

# 创建颜色
Age.col <- c("grey80",black)
Gender.col <- c(yellow, orange)
Stage.col <- alpha(blue, c(0.4, 0.6, 0.8, 1))
T.col <- alpha(green, c(0.4,0.6, 0.8, 1))
N.col <- alpha(cherry, c(0.4,0.6, 0.8, 1))
Smoking.col <- alpha(red, c(0.4,0.6, 0.8, 1))

# 硬核base plot一块一块画，当然也可以把其中的pie chart提取出来后期AI或者PPT拼接也是比较方便的
pdf("pieTable.pdf",width = 7, height = 4)
showLayout <- F # 默认不在最终pdf的首页显示layout结构，不过建议初次绘制的时候改为TRUE看一下，方便理解

# 设置画面布局，相同数字代表同一区块，数字越多代表该区块所占面积越大（一共25个区域）
layout(matrix(c( 1, 1, 1,  2, 2, 2,  3, 3, 3,  4, 4, 4,  5, 5, 5,  6, 6, 6,7, 7, 7,
                 8, 8, 8,  9, 9, 9, 10,10,10, 11,11,11, 12,12,12, 13,13,13,14,14,14,
                 8, 8, 8,  9, 9, 9, 10,10,10, 11,11,11, 12,12,12, 13,13,13,14,14,14,
                 15,15,15, 16,16,16, 17,17,17, 18,18,18,19,19,19, 20,20,20, 21,21,21,
                 15,15,15, 16,16,16, 17,17,17, 18,18,18,19,19,19, 20,20,20, 21,21,21,
                 22,22,22, 23,23,23, 24,24,24,25,25,25, 26, 26,26,27,27,27,28,28,28,
                 29,29,29, 29,29,29, 29,29,29, 29,29,29, 29,29,29, 29,29,29, 29,29,29),
              byrow = T,nrow = 7))

if(showLayout) {
  layout.show(n = 29) # 直观展示画布分布
}

#-------------------------#
# 画布区域1-6：绘制图抬头 #
#-------------------------#

par(bty="n", mgp = c(0,0,0), mar = c(0,0,0,0), lwd = 2) # 基础参数，各边界距离为0
plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "LUAD",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Age",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Gender",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Stage",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "T",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "N",cex = 1.5, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Smoking",cex = 1.5, col = "white") # 显示图标题

#--------------------------------------#
# 画布区域7-12：绘制High组抬头和扇形图 #
#--------------------------------------#

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "High\n(n = 206)",cex = 1.5, col = "white") # 显示图标题

# High group
pie(pie.high$Age$Pct, 
    col = Age.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$Gender$Pct, 
    col = Gender.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$Stage$Pct, 
    col = Stage.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$T$Pct, 
    col = T.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$N$Pct, 
    col = N.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$Smoking$Pct, 
    col = Smoking.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

abline(v = par("usr")[2], col = "black") # 右侧封上黑线

#--------------------------------------#
# 画布区域13-18：绘制Low组抬头和扇形图 #
#--------------------------------------#

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Low\n(n = 223)",cex = 1.5, col = "white") # 显示图标题

# Low group
pie(pie.low$Age$Pct, 
    col = Age.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$Gender$Pct, 
    col = Gender.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$Stage$Pct, 
    col = Stage.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$T$Pct, 
    col = T.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$N$Pct, 
    col = N.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$Smoking$Pct, 
    col = Smoking.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

abline(v = par("usr")[2], col = "black") # 右侧封上黑线

#--------------------------------#
# 画布区域19-24：绘制空抬头和p值 #
#--------------------------------#

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#333333") # 背景涂黑

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["Age"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black") # 底部封上黑线

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["Gender"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n",# 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["Stage"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["T"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["N"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black") # 底部封上黑线

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["Smoking"]),cex = 1.3, col = "#333333") # 显示图标题
abline(h = par("usr")[4], col = "black") # 底部封上黑线

abline(v = par("usr")[2], col = "black") # 右侧封上黑线

#----------------------#
# 画布区域25：绘制图例 #
#----------------------#

plot(0,0,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
legend("topleft",
       legend = c("<=65",">65",
                  "Female","Male",
                  "I","II","III","IV",
                  "T1","T2","T3","T4",
                  "N0","N1","N2","N3",
                  "False","True"),
       fill = c(Age.col,
                Gender.col,
                Stage.col,
                T.col,
                N.col,
                Smoking.col),
       border = NA, # 图例颜色没有边框
       bty = "n", # 图例没有边框
       cex = 1,
       #box.lwd = 3,
       x.intersp = 0.05,
       y.intersp = 1,
       text.width = 0.063, # 图例的间隔
       horiz = T) # 图例水平放置
abline(h = par("usr")[1], col = "black") # 底部封上黑线
abline(h = par("usr")[3], col = "black") # 底部封上黑线
abline(v = par("usr")[2], col = "black") # 右侧封上黑线
abline(h = par("usr")[4], col = "black") # 底部封上黑线
# 关闭图像句柄
invisible(dev.off())



