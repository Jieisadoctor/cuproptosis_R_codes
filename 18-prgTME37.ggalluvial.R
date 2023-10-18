#install.packages("ggplot2")
#install.packages("ggalluvial")


#???ð?
library(ggalluvial)
library(ggplot2)
library(dplyr)

riskFile="tcgaClinical.txt"            #?????ļ?
setwd("C:\\Users\\86188\\Desktop\\Bladder Endoplasmic Reticulum Stress lncRNA\\18.Sankey")     #???ù???Ŀ¼

#??ȡ?????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$N <- as.factor(risk$N)
risk$riskScore <- factor(risk$riskScore, levels = c("Low","High"))

#׼??ɣ??ͼ?????ļ?
corLodes=to_lodes_form(risk, axes = 1:ncol(risk), id = "Cohort")

#?õ??????ļ?
pdf(file="ggalluvial.pdf", width=7, height=5)
mycol=rep(c("#0066FF","#FF9900","#FF0000","#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  	 scale_x_discrete(expand = c(0, 0)) +  
  	 #??aes.flow??????????ɫ??forward˵????ɫ??ǰ??????״ͼһ?£?backward˵???ͺ???????״ͼһ?¡?
  	 geom_flow(width = 2/10,aes.flow = "forward") + 
	 geom_stratum(alpha = .9,width = 2/10) +
	 scale_fill_manual(values = mycol) +
	 #size=3??????????С
	 geom_text(stat = "stratum", size = 3,color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #ȥ????????
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = FALSE)                            
dev.off()
