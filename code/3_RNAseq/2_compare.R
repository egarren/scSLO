rm(list=ls())
library(biomaRt)
library(EnhancedVolcano)
library(cowplot)
library(dplyr)
library(ggpubr)
library(Seurat)
library(data.table)
library(VennDiagram)
library(Seurat)
library(biomaRt)
library(scales)
my.ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}


scRNA<-pdl1.DE
scRNA$gene<-rownames(scRNA)
RNA<-res_pdl1
RNA<- RNA[order(RNA$pvalue),]
RNA<-data.frame(gene=rownames(RNA),log2FC=RNA$log2FoldChange,p_val_adj=RNA$padj)
RNA<-left_join(RNA,unique(tx2gene[,c("gene_id","symbol")]),by=c("gene"="gene_id"),keep=F)
RNA<-RNA[!is.na(RNA$p_val_adj),c("log2FC","p_val_adj","symbol")]
RNA$avg_log2FC = -RNA$log2FC
hgnc_list<-getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = RNA$symbol , 
                  mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="dec2021.archive.ensembl.org"), attributesL = c("hgnc_symbol"), 
                  martL = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="dec2021.archive.ensembl.org"), uniqueRows=T)
m <- match(RNA$symbol, hgnc_list$MGI.symbol)
RNA<-cbind(RNA,hgnc_list[m,])
RNA$gene<-RNA$HGNC.symbol

#FC vs FC
df<-left_join(scRNA[,c("avg_log2FC","gene")],RNA[,c("avg_log2FC","gene")],by="gene",keep=F)
colnames(df)<-c("scRNA_FC","gene","bulkRNA_FC")
rownames(df)<-make.names(df$gene,unique=T)
df$diff<-abs(df$scRNA_FC - df$bulkRNA_FC)
df<-df[order(-df$diff),]
gene.list<-c(head(rownames(df)),"CD274")
df <- df %>% mutate(label = ifelse(gene %in% gene.list, "sig", "notsig"))
p2 <- ggplot(df, aes(bulkRNA_FC,scRNA_FC)) + geom_point(aes(colour=label))+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) +
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
        axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
  geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))+#geom_abline(intercept = 0, slope = 1)+
  scale_color_manual(values = c("sig" = hue_pal()(6)[5], "notsig" = "grey70"))
if(length(gene.list)>0){p2 <- LabelPoints(plot = p2, points = gene.list, repel = TRUE,xnudge=0,ynudge=0,size=2.5,segment.size=0.1)}
p2
ggsave2("FCvsFC.png",width=3, height=3,device="png")



#mouse vs human
m564<-data.frame(gene=DE.res_564$symbol,m564.log2FC=DE.res_564$log2FC)
wt<-data.frame(gene=DE.res_wt$symbol,wt.log2FC=DE.res_wt$log2FC)

#venn diagram
DE.list<-list()
for(i in c("DE.res_wt","DE.res_564")){
  df<-get(i)
  DE<-df[abs(df$log2FC)>2&df$p_val_adj<0.05,]
  DE.list[[i]]<-DE$symbol[!is.na(DE$symbol)]
}
names(DE.list)<-c("WT","564Igi")
venn.diagram(
  x = DE.list,category.names = names(DE.list),filename = "output/564vsWT.png",
  imagetype="png",height=700,width=700,margin=0.1,
  lwd=1,col=c("darkmagenta", "chartreuse4"),
  fill=c(alpha("darkmagenta",0.3), alpha('chartreuse4',0.3)),
  # lty="blank",fill=sample(brewer.pal(9, "Set1"),size=length(DE.list)),
  cex=0.3,fontface="bold",fontfamily="sans",
  cat.cex=0.3,cat.fontface="bold",cat.default.pos="outer",cat.fontfamily="sans",
  cat.col = c("darkmagenta", "chartreuse4")#,cat.dist = c(0.1, 0.1)#,cat.pos=c(0,180,0)
)
intersect(DE.list[[1]],DE.list[[2]]) 

#correlation
comp<-merge(wt,m564,by="gene")
rownames(comp)<-make.unique(as.character(comp$gene))
colnames(comp)<-c("gene","wt","m564")
comp<-comp[!is.na(comp$wt),]
comp<-comp[!is.na(comp$m564),]
comp<-comp[!is.na(comp$gene),]
comp$sig<-"unsig"
comp$sig[abs(comp$wt)>3&abs(comp$m564)>10]<-"co.DE"
comp$sig[abs(comp$wt)>3&abs(comp$m564)<=10]<-"wt.DE"
comp$sig[abs(comp$wt)<=3&abs(comp$m564)>10]<-"564.DE"
comp$sum<-comp$wt+comp$m564
comp$diff<-abs(comp$wt-comp$m564)
up<-head(comp$gene[order(comp$sum)],n=15)
down<-tail(comp$gene[order(comp$sum)],n=15)
diff<-head(comp$gene[order(comp$diff)],n=30)
gene.list<-unique(c(up,down,diff))
comp$sig <- factor(comp$sig, levels = c("unsig","co.DE","wt.DE","564.DE"))
p2 <- ggplot(comp, aes(wt,m564)) + geom_point(aes(colour=sig),size=0.5)+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) + 
  labs(x=bquote(~Log[2]~ (frac("Cxcl13.PDL1+SRBC","WT+SRBC"))),y=bquote(~Log[2]~ (frac("564homo.Cxcl13.PDL1","564homo"))))+theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
        axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
  geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))+#geom_abline(intercept = 0, slope = 1)+
  scale_color_manual(values = c("grey","black","red","orange"))
if(length(gene.list)>0 & length(gene.list)<=60){p2 <- LabelPoints(plot = p2, points = gene.list, repel = TRUE,xnudge=0,ynudge=0,size=2.5,segment.size=0.1)}
p2
ggsave2("output/wt.m564.scatter.png",width=5, height=5,device="png")
write.csv(comp,file="output/comp.csv")





