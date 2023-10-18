#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)       #引用包
setwd("C:\\Users\\86188\\Desktop\\LUAD-Cu\\03.maftool")      #设置工作目录

#读取突变基因文件
geneRT=read.table("gene.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneRT)

#绘制瀑布图
maf=read.maf(maf="input.maf")
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
pdf(file="oncoplot.pdf", width=9, height=6)
oncoplot(maf=maf, genes=gene, colors = vc_cols, draw_titv=T)
dev.off()

