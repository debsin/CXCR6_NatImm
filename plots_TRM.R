


getSignature<-function(fpkm,gene.list, norm = F){
  mat_interest = fpkm[intersect(gene.list, rownames(fpkm)) , , drop = F] 
  if(norm== T){
    mat_interest =log2(mat_interest[complete.cases(mat_interest),, drop = F]+1)
  }
  mat_interest = mat_interest[abs(rowSums(mat_interest))>0, ,drop = F]
  signature = apply(scale(t(mat_interest), center = F),1, mean)
  return(signature)
}

# TRM/ CXCR6 path
gene.path = "C:/Projects/external/jem_20190249_tables1/"

## tumor TRM versus. tumor non-TRM
trm_table=readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS4.xlsx"),skip = 2,col_names = T,sheet = 1)
head(trm_table)
TRMtumor = trm_table[,c("Gene ID","log2 fold change", "P adj")]
TRMtumor = TRMtumor[TRMtumor$`log2 fold change`>0 &
                      TRMtumor$`P adj`<0.05, ]$`Gene ID`

## Lung TRM vs non TRM
trm_table=readxl::read_xlsx(file.path(gene.path,"JEM_20190249_TableS2.xlsx"),skip = 2,col_names = T,sheet = 1)
TRMlung = trm_table[,c("Gene ID","log2 fold change", "P adj")]
TRMlung = TRMlung[TRMlung$`log2 fold change`>0 &
                    TRMlung$`P adj`<0.05, ]$`Gene ID`


TRMlung
TRMtumor

# Plot Literature Data
TRMlung

plot.meta =  data.frame(cbind(`Literature TRM Lung` = getSignature(assay(ntd), TRMlung),
                              `Literature TRM Lung Tumour` = getSignature(assay(ntd),TRMtumor)))

plot.meta = cbind(plot.meta, population = dds$population, patient = dds$patient)
plot.meta$patient = as.factor(dds$patient)
levels(plot.meta$population) = c("CD8+CD103-", "CD8+CD103+")


plot.data.lit = reshape::melt(plot.meta)
head(plot.data.lit)

pdf("TRM/plots/LiteratureTRM.pdf", height = 5, width = 7)
ggplot(plot.data.lit, aes(x = population, y = value, fill = population))+facet_wrap(~variable,ncol = 4)+
  geom_boxplot(width = 0.4, outlier.shape = NA)+
  geom_jitter(width = 0.07)+ylab("Expression of TRM Signature")+
  scale_fill_ipsum()+
  theme_classic()+theme(legend.position = "none", axis.title.x = element_blank())+
  ggtitle("Human Tumour Patients")+
  stat_compare_means(comparisons = list( c("CD8+CD103-", "CD8+CD103+")), 
                     label.y = 1.2,  label = "p.format", 
                     method = "wilcox.test",  size = 3.5)
dev.off()



mice.tmem.v.rest = rownames(tmem.v.rest)[tmem.v.rest$avg_log2FC>2 & tmem.v.rest$p_val_adj<0.05] 
mice.tmem.v.resthomo = unique(mus2humanGENE(mice.tmem.v.rest))


plot.meta =  data.frame(cbind(`Mice TRM` = getSignature(assay(ntd), mice.tmem.v.resthomo),
                              population = dds$population, patient = dds$patient))
plot.meta$patient = as.factor(dds$patient)
plot.meta$population = as.factor(dds$population)
levels(plot.meta$population) = c("CD8+CD103-", "CD8+CD103+")


plot.data.mice = reshape::melt(plot.meta)
head(plot.data.mice)

pdf("TRM/plots/MiceTRM.pdf", height = 5, width = 5)
ggplot(plot.data.mice, aes(x = population, y = value, fill = population))+
  geom_boxplot(width = 0.4, outlier.shape = NA)+
  geom_jitter(width = 0.07)+ylab("Expression of Signature")+
  scale_fill_ipsum()+
  theme_classic()+theme(legend.position = "none", axis.title.x = element_blank())+
  ggtitle("Human Tumour Patients", subtitle = "TREM geneset derived from Mice data")+
  stat_compare_means(comparisons = list( c("CD8+CD103-", "CD8+CD103+")), 
                     label.y = 1.15,  label = "p.format", 
                     method = "wilcox.test",  size = 3.5)
dev.off()



df <-data.frame(colData(dds)[,c("patient", "population"), drop= F])
order.idx = order(df$population)
df = df[order.idx, ,drop = F]
select = which(rownames(assay(ntd)) %in% mice.tmem.v.resthomo)

annoCol<-list(population=c(CD8posCD103neg ="blue", CD8posCD103pos ="red"))


pdf("TRM/plots/heatmapHuman.pdf", height = 6, width = 6)

pheatmap(assay(ntd)[select,order.idx], 
         # clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         cluster_rows=T, cluster_cols = T,scale = "row",
         annotation_col=df, annotation_colors = annoCol)
dev.off()


pdf("TRM/plots/umapMiceTRM.pdf", height = 5.5, width = 6)
FeaturePlot(mice, features = "mice.tmem.v.rest1", min.cutoff = 0.1)+ ggtitle("Tmem (vs Rest) Signature", subtitle = "from mice")
dev.off()

pdf("TRM/plots/umapMiceCT.pdf", height = 5.5, width = 6)
DimPlot(mice, group.by = "annotation", label = T, order = T, label.box = T)+ggtitle("T mem Annotaion")+NoLegend()
dev.off()

pdf("TRM/plots/umapHumanTRM.pdf", height = 5.5, width = 6)
FeaturePlot(mice, features = "TRM_human_expt1", min.cutoff = 0.05, order = T)+ggtitle("TRM Human Homolog", subtitle = "from HUMAN CD8+ tumour experiment")
dev.off()


# convert mice TRM and project on human
human.TRM_homoT0 = unique(mus2humanGENE(rownames(trm.mice.t0)[trm.mice.t0$avg_log2FC>0&trm.mice.t0$p_val_adj<0.05]))
human.TRM_homoT1 = unique(mus2humanGENE(rownames(trm.mice.t1)[trm.mice.t1$avg_log2FC>0&trm.mice.t1$p_val_adj<0.05]))
human.TRM_homoT3 =unique(mus2humanGENE( rownames(trm.mice.t3)[trm.mice.t3$avg_log2FC>0&trm.mice.t3$p_val_adj<0.05]))
human.TRM_homoT7 = unique(mus2humanGENE(rownames(trm.mice.t7)[trm.mice.t7$avg_log2FC>0&trm.mice.t7$p_val_adj<0.05]))

plot.meta =  data.frame(cbind(miceTRMT0 = getSignature(assay(ntd), human.TRM_homoT0),
                              miceTRMT1 = getSignature(assay(ntd), human.TRM_homoT1),
                              miceTRMT3 = getSignature(assay(ntd), human.TRM_homoT3),
                              miceTRMT7 = getSignature(assay(ntd),  human.TRM_homoT7)))

plot.meta = cbind(plot.meta, population = dds$population, patient = dds$patient)
plot.meta$patient = as.factor(dds$patient)
plot.meta$population = as.factor(dds$population)
levels(plot.meta$population) = c("CD8+CD103-", "CD8+CD103+")


plot.data.tmem.tx = reshape::melt(plot.meta)
head(plot.data)

vcn_panel = ggplot(plot.data.tmem.tx, aes(x = population, y = value, fill = population))+facet_wrap(~variable,ncol = 4)+
  geom_boxplot(width = 0.4, outlier.shape = NA)+
  geom_jitter(width = 0.07)+ylab("Expression of Signature")+
  scale_fill_ipsum()+
  theme_classic()+theme(legend.position = "none", axis.title.x = element_blank())+
  ggtitle("Human Tumour patients", subtitle = "Signature genesets derived from mice data (each.Tx-vs-rest.Tx) within Tmem")+
  stat_compare_means(label = "p.format", method = "wilcox.test", vjust = 1, size = 3.5)

p0 = plotVolCano(trm.mice.t0)+ggtitle("T0 vs rest Tx within Tmem")
p1 = plotVolCano(trm.mice.t1)+ggtitle("T1 vs rest Tx within Tmem")
p3 = plotVolCano(trm.mice.t3)+ggtitle("T3 vs rest Tx within Tmem")
p7 = plotVolCano(trm.mice.t7)+ggtitle("T7 vs rest Tx within Tmem")


# Pivot: Mice one.Tx-vs-rest.Tx points within Tmem
pdf("TRM/plots/Tx Mice TRM.pdf", height = 6, width = 15)

gridExtra::grid.arrange(p0,p1,p3,p7, ncol = 4)
vcn_panel

dev.off()




## Pivot: Mice Tmem vs Rest within each Tx
homolog.T0.TRM = unique(mus2humanGENE(rownames(t0.mice.trm)[t0.mice.trm$avg_log2FC>0&t0.mice.trm$p_val_adj<0.05]))
homolog.T1.TRM = unique(mus2humanGENE(rownames(t1.mice.trm)[t1.mice.trm$avg_log2FC>0&t1.mice.trm$p_val_adj<0.05]))
homolog.T3.TRM =unique(mus2humanGENE( rownames(t3.mice.trm)[t3.mice.trm$avg_log2FC>0&t3.mice.trm$p_val_adj<0.05]))
homolog.T7.TRM = unique(mus2humanGENE(rownames(t7.mice.trm)[t7.mice.trm$avg_log2FC>0&t7.mice.trm$p_val_adj<0.05]))

plot.meta =  data.frame(cbind(miceTRMT0 = getSignature(assay(ntd), homolog.T0.TRM),
                              miceTRMT1 = getSignature(assay(ntd), homolog.T1.TRM),
                              miceTRMT3 = getSignature(assay(ntd), homolog.T3.TRM),
                              miceTRMT7 = getSignature(assay(ntd), homolog.T3.TRM)))

plot.meta = cbind(plot.meta, population = dds$population, patient = dds$patient)
plot.meta$patient = as.factor(dds$patient)
plot.meta$population = as.factor(dds$population)
levels(plot.meta$population) = c("CD8+CD103-", "CD8+CD103+")


plot.data.tmem.tx = reshape::melt(plot.meta)
head(plot.data)

vcn_panel = ggplot(plot.data.tmem.tx, aes(x = population, y = value, fill = population))+facet_wrap(~variable,ncol = 4)+
  geom_boxplot(width = 0.4, outlier.shape = NA)+
  geom_jitter(width = 0.07)+ylab("Expression of Signature")+
  scale_fill_ipsum()+
  theme_classic()+theme(legend.position = "none", axis.title.x = element_blank())+
  ggtitle("Human Tumour patients", subtitle = "Signature genesets derived from mice data (Tmem-vs-rest T) within each Tx")+
  stat_compare_means(label = "p.format", method = "wilcox.test", vjust = 1, size = 3.5)

p0 = plotVolCano(t0.mice.trm)+ggtitle("Tmem vx rest T within T0")
p1 = plotVolCano(t1.mice.trm)+ggtitle("Tmem vx rest T within T1")
p3 = plotVolCano(t3.mice.trm)+ggtitle("Tmem vx rest T within T3")
p7 = plotVolCano(t7.mice.trm)+ggtitle("Tmem vx rest T within T7")


# Pivot: Mice one.Tx-vs-rest.Tx points within Tmem
pdf("TRM/plots/Tmem Mice Tx.pdf", height = 6, width = 15)

gridExtra::grid.arrange(p0,p1,p3,p7, ncol = 4)
vcn_panel

dev.off()
