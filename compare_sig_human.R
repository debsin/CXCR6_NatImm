library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
library(pheatmap)
library(ggpubr)

human_t = read.csv("C:/Projects/CXCR6/TRM/HUMAN TRM TUMOR/P26_MALE_Expression_table.tsv", sep = "\t", header = T, row.names = 1)
human_t.meta = as.data.frame(readxl::read_xlsx("C:/Projects/CXCR6/TRM/HUMAN TRM TUMOR/Sample_ID_RNAseq TRM.xlsx"))

rownames(human_t.meta) = human_t.meta$`Sample ID`
human_t.meta = human_t.meta[, -1]


mart_export <- read.delim("C:/Projects/CXCR6/TRM/HUMAN TRM TUMOR/mart_export.txt")

mart_export = mart_export[, c("Gene.stable.ID", "Transcript.stable.ID","Transcript.stable.ID.version", "Gene.name")]
mart_export = mart_export[!duplicated(mart_export), ]
rownames(mart_export) = mart_export$Transcript.stable.ID.version
mart_export = mart_export[mart_export$Gene.name != "",]
mart_export = mart_export[complete.cases(mart_export),]

# mart_export = mart_export[!duplicated(mart_export$Gene.name), ]
a = data.frame(table(mart_export$Transcript.stable.ID.version))
a = a[order(a$Freq, decreasing = T),]
head(a)

t(mart_export[which(mart_export$Transcript.stable.ID.version == "ENST00000636378.2"), ])


matching_rows = intersect(rownames(human_t) , rownames(mart_export))
mart_export = mart_export[matching_rows, ]
human_match = human_t[matching_rows, ]
# rownames(human_match) = mart_export[matching_rows,]$Gene.name

which(mart_export$Gene.name=="A1BG")
human_match[mart_export$Transcript.stable.ID.version[which(mart_export$Gene.name=="A1BG")], ]

human_match$Gene = mart_export[rownames(human_match),]$Gene.name


human_raw = aggregate(.~Gene,data=human_match,FUN=sum)
rownames(human_raw) = human_raw$Gene
human_raw = round(human_raw[,-1])


colnames(human_raw)
coldata = human_t.meta
coldata$population= gsub("[+]", "pos", coldata$population)
coldata$population = gsub("[-]", "neg", coldata$population)

all(rownames(coldata) %in% colnames(human_raw))
# human_raw <- human_raw[, rownames(coldata)]

all(rownames(coldata) == colnames(human_raw))

apply(human_raw, 2, class)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = human_raw,
                              colData = coldata,
                              design = ~patient+population)
dds


keep <- rowSums(counts(dds)) >= 10

table(keep)

dds <- dds[keep,]


dds$population <- relevel(dds$population, ref = "CD8posCD103neg")
dds <- DESeq(dds)
res <- results(dds)
res
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

# res05 <- results(dds, alpha=0.05)
# summary(res05)


ntd <- normTransform(dds)

vsd <- vst(dds, blind=FALSE)



library(ggrepel)

plotVolCano<-function(genes){
  # genes$padj = genes$padj
  # genes$log2FoldChange = genes$log2FoldChange
  # genes$pvalue  = genes$p_val
  genes$Gene = rownames(genes)
  genes[which(genes[,"pvalue"] == 0), "pvalue"] <- .Machine$double.xmin
  
  genes$Significant <- ifelse(genes$padj < 0.05 & abs(genes$log2FoldChange) > 0, "FDR < 0.05", "Not Sig")
  
  ggplot(genes, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = Significant)) +
    scale_color_manual(values = c("red", "grey")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_text_repel(
      data = subset(genes, genes$Significant != "Not Sig"),
      aes(label = Gene),
      size = 3.1,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )
}

plotVolCano(as.data.frame(res))+ggtitle("CD8+CD103+ vs CD8+CD103-")


## Define TRM signature in Tumour
TRMup = rownames(res)[which((res$log2FoldChange>0 & res$padj<0.05) == T)]
TRMdn = rownames(res)[which((res$log2FoldChange<0 & res$padj<0.05) == T)]

plotCounts(dds, gene="ITGAE", intgroup="population", pch = 16)


getSignature<-function(fpkm,gene.list, norm = F){
  mat_interest = fpkm[intersect(gene.list, rownames(fpkm)) , , drop = F] 
  if(norm== T){
    mat_interest =log2(mat_interest[complete.cases(mat_interest),, drop = F]+1)
  }
  mat_interest = mat_interest[abs(rowSums(mat_interest))>0, ,drop = F]
  signature = apply(scale(t(mat_interest), center = F),1, mean)
  return(signature)
}

windowsFonts("Arial" = windowsFont("Arial"))

test.df = data.frame(TRM_raw = getSignature(assay(dds), TRMup))
test.df = cbind(test.df, TRM_norm = getSignature(assay(ntd), TRMup))
test.df = cbind(test.df, population = dds$population)

test.df= melt(test.df)
head(test.df)

ggplot(data=test.df, aes(x=value, group=population, fill=population)) +
  geom_density(adjust=1.5, alpha=.8) + facet_wrap(~variable)+scale_fill_ipsum()+
  theme_ipsum()


plot.meta = data.frame(TRMup = getSignature(assay(ntd), TRMup))
# plot.meta$TRMdn = getSignature(assay(ntd), TRMdn)
plot.meta$TRMlung = getSignature(assay(ntd), TRMlung)
plot.meta$TRMtumor = getSignature(assay(ntd), TRMtumor)

plot.meta = cbind(plot.meta, population = dds$population)

plot.meta$patient = as.factor(dds$patient)
plot.meta$population = as.factor(plot.meta$population)

plot.data = reshape::melt(plot.meta)
head(plot.data)

g1 = ggplot(data=plot.data, aes(x=value, group=population, fill=population)) +
  geom_density(adjust=1.5, alpha=.9) + facet_wrap(~variable)+scale_fill_ipsum()+
  theme_ipsum()+theme(legend.position = "bottom")+ggtitle("Human Tumour Cells")

g2 = ggplot(plot.data, aes(x = population, y = value, fill = population))+facet_wrap(~variable)+
  geom_boxplot(width = 0.4, outlier.shape = NA)+
  geom_jitter(width = 0.05)+ylab("TRM Signature")+
  theme_ipsum()+theme(legend.position = "bottom", axis.text.x=element_blank())+
  scale_fill_ipsum()+
  stat_compare_means(label = "p.format", method = "wilcox.test", vjust = 0.95, size = 3.5)


gridExtra::grid.arrange(g1,g2, ncol = 1)




# find  human homologs in mouse

mice.TRMup_homo = unique(human2musGENE(TRMup))
mice.TRMdn_homo = unique(human2musGENE(TRMdn))


mice = AddModuleScore(mice, features = list("TRM_human_expt" = mice.TRMup_homo), 
                      name = "TRM_human_expt", seed = 10,assay = "RNA")




dittoPlot(mice, var = "TRM_human_skin1", group.by = "reducedCellType")+geom_boxplot(fill="white", width = 0.5)+
  ggtitle("TRM Signature from HUMAN Dermis(JEM paper)")+NoLegend()+
  dittoPlot(mice.alexis, var = "TRM_human_skin1", group.by = "celltype")+geom_boxplot(fill="white", width = 0.5)+
  ggtitle("Reduced Annotation")+NoLegend()


mice.alexis  = subset(mice, cells = grep("Alexis", mice$celltype))

dittoDotPlot(mice.alexis, vars = intersect(mice.TRMup_homo, rownames(mice.alexis)), group.by = "reducedCellType" )+
  xlab("Found homologs")+coord_flip()

length(TRMup) - length(intersect(mice.TRMup_homo, rownames(mice.alexis)))


# mice$tmem_group = mice$celltype
# mice$tmem_group[grep("Tmem", mice$tmem_group, invert = T)] = "others"



f1 = FeaturePlot(mice, features = "TMEM1", min.cutoff = 0.1, order = T)+ggtitle("T.mem signature common Tx", subtitle = "from mice")
f2 = FeaturePlot(mice, features = "mice.tmem.v.rest1", min.cutoff = 0.1)+ ggtitle("Tmem vs Rest Signature", subtitle = "from mice")

f3 = FeaturePlot(mice, features = "TRM_human_Lung1", min.cutoff = 0.05, order = T)+ggtitle("TRM Human Lung", subtitle ="from literature")
f4 = FeaturePlot(mice, features = "TRM_human_Tumor1", min.cutoff = 0.05, order = T)+ggtitle("TRM  Human Lung Tumor", subtitle ="from literature")

f5 = FeaturePlot(mice, features = "TRM_human_expt1", min.cutoff = 0.05, order = T)+ggtitle("TRM Human Homolog", subtitle = "from HUMAN CD8+ tumour experiment")
d1 = DimPlot(mice, group.by = "celltype", label = T, label.size = 2, label.box = T)+NoLegend()

vn1 = dittoPlot(mice, var = "TMEM1", group.by = "celltype" )+
  geom_boxplot(fill = "white", width = 0.5)+ggtitle("T.mem signature common Tx", subtitle = "from mice")+
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref = "Alexis.Tmem")+NoLegend()

vn2 = dittoPlot(mice, var = "mice.tmem.v.rest1", group.by = "celltype" )+
  geom_boxplot(fill = "white", width = 0.5)+ggtitle("Tmem vs Rest Signature", subtitle = "from mice")+
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref = "Alexis.Tmem")+NoLegend()

vn3 = dittoPlot(mice, var = "TRM_human_Lung1", group.by = "celltype" )+
  geom_boxplot(fill = "white", width = 0.5)+ggtitle("TRM Human Lung", subtitle ="from literature")+
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref = "Alexis.Tmem")+NoLegend()

vn4 = dittoPlot(mice, var = "TRM_human_Tumor1", group.by = "celltype" )+
  geom_boxplot(fill = "white", width = 0.5)+ggtitle("TRM  Human Lung Tumor", subtitle ="from literature")+
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref = "Alexis.Tmem")+NoLegend()

vn5 = dittoPlot(mice, var = "TRM_human_expt1", group.by = "celltype" )+
  geom_boxplot(fill = "white", width = 0.5)+ggtitle("TRM Human Homolog", subtitle = "from HUMAN CD8+ tumour experiment")+
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref = "Alexis.Tmem")+NoLegend()


gridExtra::grid.arrange(f1,f2, f3,f4,f5,
                        vn1,vn2, vn3,vn4,vn5,ncol = 5)






# convert mice TRM and project on human
human.TRM_homo = unique(mus2humanGENE(mice.tmem.v.rest))
human.TRM_homoT0 = unique(mus2humanGENE(rownames(trm.mice.t0)[trm.mice.t0$avg_log2FC>0&trm.mice.t0$p_val_adj<0.05]))
human.TRM_homoT1 = unique(mus2humanGENE(rownames(trm.mice.t1)[trm.mice.t1$avg_log2FC>0&trm.mice.t1$p_val_adj<0.05]))
human.TRM_homoT3 =unique(mus2humanGENE( rownames(trm.mice.t3)[trm.mice.t3$avg_log2FC>0&trm.mice.t3$p_val_adj<0.05]))
human.TRM_homoT7 = unique(mus2humanGENE(rownames(trm.mice.t7)[trm.mice.t7$avg_log2FC>0&trm.mice.t7$p_val_adj<0.05]))

miceTRMhomo = intersect(rownames(dds), mus2humanGENE(micecl17homo))


plotCounts(dds, gene = , intgroup="population", pch = 16)


plot.meta =  data.frame(cbind(TRMup = getSignature(assay(ntd), TRMup),
                              LitTRMlung = getSignature(assay(ntd), TRMlung),
                              LitTRMtumor = getSignature(assay(ntd),TRMtumor),
                              miceTRM = getSignature(assay(ntd), human.TRM_homo),
                              # micecl17homo = getSignature(assay(ntd), micecl17homo),
                              miceTRMT0 = getSignature(assay(ntd), human.TRM_homoT0),
                              miceTRMT1 = getSignature(assay(ntd), human.TRM_homoT1),
                              miceTRMT3 = getSignature(assay(ntd), human.TRM_homoT3),
                              miceTRMT7 = getSignature(assay(ntd),  human.TRM_homoT7)))

plot.meta = cbind(plot.meta, population = dds$population, patient = dds$patient)
plot.meta$patient = as.factor(dds$patient)
plot.meta$population = as.factor(plot.meta$population)

plot.data.tmem.tx = reshape::melt(plot.meta)
head(plot.data)

ggplot(plot.data, aes(x = population, y = value, fill = population))+facet_wrap(~variable,ncol = 4)+
  geom_boxplot(width = 0.4, outlier.shape = NA)+
  geom_jitter(width = 0.05)+ylab("TRM Signature from Mice")+
  theme_classic()+theme(legend.position = "none")+ggtitle("Human Tumour patients")+
  stat_compare_means(label = "p.format", method = "wilcox.test", vjust = 1, size = 3.5)



plotPCA(vsd, intgroup=c("population", "patient"))


pcaData <- plotPCA(vsd, intgroup=c("population", "patient"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcp1 = ggplot(pcaData, aes(PC1, PC2, color=population, label=patient)) +
  geom_point(size=3) + geom_text_repel(size = 3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+theme_classic()+theme(legend.position = "bottom")+ggtitle("all features")


top_mice.trm
miceTRMhomo = intersect(rownames(dds), mus2humanGENE(top_mice.trm))

select = which(rownames(assay(vsd)) %in% miceTRMhomo)

pcaData <- plotPCA(vsd[select,], intgroup=c("population", "patient"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcp2 = ggplot(pcaData, aes(PC1, PC2, color=population, label=patient)) +
  geom_point(size=3) + geom_text_repel(size = 3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()+
  theme_classic()+theme(legend.position = "bottom")+ggtitle("Mice TRM")


df <-data.frame(colData(dds)[,c("patient", "population"), drop= F])
order.idx = order(df$population)
df = df[order.idx, ,drop = F]
select = which(rownames(assay(vsd)) %in% miceTRMhomo)

annoCol<-list(population=c(CD8posCD103neg ="blue", CD8posCD103pos ="red"))
pheatmap(assay(ntd)[select,order.idx], 
         cluster_rows=T, cluster_cols = T,scale = "row",
         annotation_col=df, annotation_colors = annoCol)



## Mice TRM vs REST
tmem.v.rest
mice.tmem.v.rest = rownames(tmem.v.rest)[tmem.v.rest$avg_log2FC>2 & tmem.v.rest$p_val_adj<0.05] 
grep("Cxcr6", mice.tmem.v.rest, value = T)
mice.tmem.v.resthomo = intersect(rownames(dds), mus2humanGENE(mice.tmem.v.rest))
mice.tmem.v.resthomo.sig = getSignature(assay(ntd), mice.tmem.v.resthomo)

feature.df = data.frame("MiceTMEMvRest" = mice.tmem.v.resthomo.sig, "population" = dds$population)
head(feature.df)
feature.df = melt(feature.df)
head(feature.df)

PT.trm.human.mice = 
  ggplot(feature.df, aes(x = population, y= value, group = population, fill = population))+ 
  geom_boxplot(width = 0.4, outlier.shape = NA)+
  geom_jitter(width = 0.05, size = 1)+
  theme_classic()+theme(legend.position = "bottom", axis.text.x=element_blank())+
  ggtitle("TRM genes from mice on Human Tumour patients")+
  stat_compare_means(label = "p.format", method = "wilcox.test", vjust = 1, size = 3.5)

PT.trm.human.mice


select = which(rownames(assay(vsd)) %in% mice.tmem.v.resthomo)

pcaData <- plotPCA(ntd[select,], intgroup=c("population", "patient"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=population, label=patient)) +
  geom_point(size=3) + geom_text_repel(size = 3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()+
  theme_classic()+theme(legend.position = "bottom")+ggtitle("Mice TRM")



df <-data.frame(colData(dds)[,c("patient", "population"), drop= F])
order.idx = order(df$population)
df = df[order.idx, ,drop = F]
select = which(rownames(assay(vsd)) %in% mice.tmem.v.resthomo)

annoCol<-list(population=c(CD8posCD103neg ="blue", CD8posCD103pos ="red"))
pheatmap(assay(ntd)[select,order.idx], 
         # clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         cluster_rows=T, cluster_cols = T,scale = "row",
         annotation_col=df, annotation_colors = annoCol)




# Mice TRM + Cluster 17 ve rest


# 
# top_mice.trm
# mice.clust17 = rownames(clust17markers)[clust17markers$avg_log2FC>2 & clust17markers$p_val_adj<0.05] 
# grep("Cxcr6", mice.clust17, value = T)
# micecl17homo = intersect(rownames(dds), mus2humanGENE(mice.clust17))
# micecl17homo.sig = getSignature(assay(ntd), micecl17homo)
# 
# feature.df = data.frame("MiceClust17" = micecl17homo.sig, "population" = dds$population)
# head(feature.df)
# feature.df = melt(feature.df)
# head(feature.df)
# 
# ggplot(feature.df, aes(x = population, y= value, group = population, fill = population))+ 
#   geom_boxplot(width = 0.4, outlier.shape = NA)+
#   geom_jitter(width = 0.05, size = 1)+
#   theme_classic()+theme(legend.position = "bottom", axis.text.x=element_blank())+
#   ggtitle("Mouse filtered TRM+cluster17 genes on Human Tumour patients")+
#   stat_compare_means(label = "p.format", method = "wilcox.test", vjust = 1, size = 3.5)
# 
# intersect(mice.tmem.v.resthomo, micecl17homo)

##
##
# miceTRMhomo = intersect(rownames(dds), mus2humanGENE(micecl17homo))
# feature.df = data.frame(t(assay(ntd)[micecl17homo, ]), "population" = dds$population)
# head(feature.df)
# feature.df = melt(feature.df)
# head(feature.df)
# 
# ggplot(feature.df, aes(x = population, y= value, group = population, fill = population))+ facet_wrap(~variable, ncol = 7)+
#   geom_boxplot(width = 0.4, outlier.shape = NA)+
#   geom_jitter(width = 0.05, size = 0.5)+
#   theme_classic()+theme(legend.position = "bottom", axis.text.x=element_blank())+
#   ggtitle("Mouse TRM genes on Human Tumour patients")+
#   stat_compare_means(label = "p.format", method = "wilcox.test", vjust = 1, size = 3.5)
# 

















micecl17homo

select = which(rownames(assay(vsd)) %in% micecl17homo)

pcaData <- plotPCA(vsd[select,], intgroup=c("population", "patient"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcp3 = ggplot(pcaData, aes(PC1, PC2, color=population, label=patient)) +
  geom_point(size=3) + geom_text_repel(size = 3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()+
  theme_classic()+theme(legend.position = "bottom")+ggtitle("Mice TRM + Cluster 17")


gridExtra::grid.arrange(pcp1,pcp2,pcp3,ncol = 3)



library("pheatmap")

miceTRMhomo = mus2humanGENE(top_mice.trm)

df <-data.frame(colData(dds)[,c("patient", "population"), drop= F])
order.idx = order(df$population)
df = df[order.idx, ,drop = F]
select = which(rownames(assay(vsd)) %in% miceTRMhomo)

annoCol<-list(population=c(CD8posCD103neg ="blue", CD8posCD103pos ="red"))
pheatmap(assay(vsd)[select,order.idx], 
         cluster_rows=T, cluster_cols = F,scale = "row",
         annotation_col=df, annotation_colors = annoCol)





df <-data.frame(colData(dds)[,c("patient", "population"), drop= F])
order.idx = order(df$population)
df = df[order.idx, ,drop = F]
select = which(rownames(assay(vsd)) %in% micecl17homo)

annoCol<-list(population=c(CD8posCD103neg ="blue", CD8posCD103pos ="red"))
pheatmap(assay(vsd)[select,order.idx], 
         cluster_rows=T, cluster_cols = F,scale = "row",
         annotation_col=df, annotation_colors = annoCol)


intersect(TRMup, TRMtumor)






























top_mice.trm
miceTRMhomo = mus2humanGENE(top_mice.trm)

box.df = data.frame(t(assay(vsd)[intersect(rownames(assay(vsd)), miceTRMhomo),]))
box.df$population = colData(vsd)[,c("population")]

box.df = reshape::melt(box.df)
head(box.df)
# pdf("genes.pdf", height = 80, width = 80)
ggplot(box.df,aes(x= population, fill = population, y = value))+
  geom_boxplot()+
  facet_wrap(~variable, ncol = 10)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",vjust = 1)+
  theme_classic()+theme(legend.position = "bottom",
                        axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank())
# 
# 
# dev.off()


#### 

# plot.meta$miceTRM = data.frame(getSignature(assay(dds), human.TRM_homo))
plot.meta =  data.frame(miceTRMupT0= getSignature(assay(ntd), human.TRM_homoT0))
plot.meta = cbind(plot.meta, population = dds$population)
plot.meta$miceTRMupT1 = getSignature(assay(ntd), human.TRM_homoT1)
plot.meta$miceTRMupT3 = getSignature(assay(ntd), human.TRM_homoT3)
plot.meta$miceTRMupT7 = getSignature(assay(ntd), human.TRM_homoT7)

plot.meta$patient = as.factor(dds$patient)
plot.meta$population = as.factor(plot.meta$population)

plot.data = reshape::melt(plot.meta)
head(plot.data)


ggplot(plot.data, aes(x = population, y = value, fill = population))+facet_wrap(~variable,ncol = 4)+
  geom_boxplot(width = 0.4, outlier.shape = NA)+
  geom_jitter(width = 0.05)+ylab("TRM Signature from Mice")+
  theme_classic()+theme(legend.position = "none")+ggtitle("Human Tumour Cells")+
  stat_compare_means(label = "p.format", method = "wilcox.test")





intersect(human.TRM_homo, TRMdn)



resultsNames(dds)
resLFC <- lfcShrink(dds, coef="population_CD8posCD103pos_vs_CD8posCD103neg", type="apeglm")


plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]


plotCounts(dds, gene=which.min(res$padj), intgroup="population")
mcols(res)$description


# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)



















## Cooks 

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

plotDispEsts(dds)

plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)


W <- res$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))

