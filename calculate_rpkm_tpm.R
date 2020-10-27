rm(list = ls())#删除目前工作目录的变量
library(xlsx)
library(readxl)
ann<- read_excel("Counts.xlsx")#读取基因文件
input<- read_excel("len.xlsx")#自己在Excel中把Data里的txt文件基因和长度共两列提取出来
library(dplyr)
merge<-left_join(ann,input,by="Gene")#根据基因那列进行合并
merge <- na.omit(merge)#删除错误值行
write.csv(merge,file = "merge.csv",sep = "\t")#写出文件

mycounts<-read.csv("merge.csv")
head(mycounts)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
head(mycounts)#最后一列Length是基因长度

#TPM计算
kb <- mycounts$Length / 1000
kb
countdata <- mycounts[,1:45]
rpk <- countdata / kb
rpk
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
write.table(tpm,file="tpm.xls",sep="\t",quote=F)

fpkm <- t(t(rpk)/colSums(countdata) * 10^6) 
head(fpkm)
write.table(fpkm,file="fpkm.xls",sep="\t",quote=F)
