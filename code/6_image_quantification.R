rm(list=ls())
library(tidyverse)
library(readxl)
library(ggpubr)
library(cowplot)
library(zoo)

metadata<-read_excel("Mouse Record.xlsx",guess_max=10000) #define path\
meta<-metadata[,c("mouse","Gender","Age","Treatment","tx2")]


df<-read.csv("test2Image.csv")
df2<-read.csv("test3Image.csv")
df<-rbind(df,df2)
df<-cbind(df$FileName_AB,df[,grepl("MeanIntensity",colnames(df))])
colnames(df)<-c("file","mem","non_mem")
samples<-as.data.frame(str_split(str_split_i(df$file,"-",5),"_",simplify=T))
df$mouse<-samples$V2
df$organ<-gsub(".tif","",samples$V3)
df$mouse<-as.numeric(df$mouse)
df<-left_join(df,meta[,c("mouse","tx2")],by="mouse")


df<-df[order(df$tx2),]
write.csv(df,"pdl1_mfi.csv")
